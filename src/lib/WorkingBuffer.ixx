export module WorkingBuffer;

import std;
import Options;

import "cdhit-common.h";

export struct WorkingBuffer
{
	std::vector<int> taap;
	std::vector<int> word_encodes;
	std::vector<int> word_encodes_backup;
	std::vector<INTs> word_encodes_no;
	std::vector<INTs> aap_list;
	std::vector<INTs> aap_begin;
	// std::vector<IndexCount>  indexCounts;
	NVector<IndexCount> lookCounts;
	NVector<uint32_t> indexMapping;
	MatrixInt64 score_mat;
	MatrixInt back_mat;
	std::vector<int> diag_score;
	std::vector<int> diag_score2;
	std::vector<int> aan_list_comp;
	std::vector<char> seqi_comp;
	int total_bytes;

	WorkingBuffer(size_t frag = 0, size_t maxlen = 0)
	{
		Set(frag, maxlen);
		seqi_comp.resize(MAX_SEQ);
	}
	void Set(size_t frag, size_t maxlen)
	{
		bool est = options.isEST;
		size_t m = MAX_UAA * MAX_UAA;
		size_t max_len = maxlen;
		size_t band = max_len * max_len;
		if (band > options.band_width)
            band = options.band_width;
		if (est)
			m = m * m;
		taap.resize(m);
		aap_list.resize(max_len);
		aap_begin.resize(m);
		// indexCounts.resize( max_len );
		word_encodes.resize(max_len);
		word_encodes_no.resize(max_len);
		word_encodes_backup.resize(max_len);
		/* each table can not contain more than MAX_TABLE_SEQ representatives or fragments! */
		if (frag > MAX_TABLE_SEQ)
			frag = MAX_TABLE_SEQ;
		lookCounts.Resize(frag + 2);
		indexMapping.Resize(frag + 2);
		diag_score.resize(MAX_DIAG);
		diag_score2.resize(MAX_DIAG);
		aan_list_comp.resize(max_len);
		total_bytes = max_len;
		total_bytes += taap.size() * sizeof(int);
		total_bytes += word_encodes.size() * sizeof(int);
		total_bytes += word_encodes_backup.size() * sizeof(int);
		total_bytes += diag_score.size() * sizeof(int);
		total_bytes += diag_score2.size() * sizeof(int);
		total_bytes += aan_list_comp.size() * sizeof(int);
		total_bytes += word_encodes_no.size() * sizeof(INTs);
		total_bytes += aap_list.size() * sizeof(INTs);
		total_bytes += aap_begin.size() * sizeof(INTs);
		total_bytes += indexMapping.Size() * sizeof(uint32_t);
		// total_bytes += indexCounts.size()*sizeof(IndexCount);
		total_bytes += lookCounts.Size() * sizeof(IndexCount);
		total_bytes += max_len * (band * sizeof(int) + sizeof(VectorInt));
		total_bytes += max_len * (band * sizeof(int) + sizeof(VectorInt64));
	}

	int EncodeWords(Sequence* seq, int NAA, bool est)
	{
		char* seqi = seq->data;
		int len = seq->size;
		// check_word_encodes
		int aan_no = len - NAA + 1;
		int i, j, i0, i1;
		int skip = 0;
		unsigned char k, k1;
		for (j = 0; j < aan_no; j++) {
			char* word = seqi + j;
			int encode = 0;
			for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) encode += word[k] * NAAN_array[k1];
			word_encodes[j] = word_encodes_backup[j] = encode;
		}

		if (est) {
			for (j = 0; j < len; j++) {
				if (seqi[j] >= 4) {                      // here N is 4
					i0 = (j - NAA + 1 > 0) ? j - NAA + 1 : 0;
					i1 = j < aan_no ? j : aan_no - 1;
					for (i = i0; i <= i1; i++) word_encodes[i] = -1;
				}
			}
			for (j = 0; j < aan_no; j++) skip += (word_encodes[j] == -1);
		}

		std::sort(word_encodes.begin(), word_encodes.begin() + aan_no);
		for (j = 0; j < aan_no; j++) word_encodes_no[j] = 1;
		for (j = aan_no - 1; j; j--) {
			if (word_encodes[j] == word_encodes[j - 1]) {
				word_encodes_no[j - 1] += word_encodes_no[j];
				word_encodes_no[j] = 0;
			}
		}
		return skip;
		// END check_word_encodes
	}
	void ComputeAAP(const char* seqi, int size)
	{
		size_t len1 = size - 1;
		for (size_t sk = 0; sk < NAA2; sk++) taap[sk] = 0;
		for (size_t j1 = 0; j1 < len1; j1++) {
			size_t c22 = seqi[j1] * NAA1 + seqi[j1 + 1];
			taap[c22]++;
		}
		for (size_t sk = 0, mm = 0; sk < NAA2; sk++) {
			aap_begin[sk] = mm; mm += taap[sk]; taap[sk] = 0;
		}
		for (size_t j1 = 0; j1 < len1; j1++) {
			size_t c22 = seqi[j1] * NAA1 + seqi[j1 + 1];
			aap_list[aap_begin[c22] + taap[c22]++] = j1;
		}
	}
	void ComputeAAP2(const char* seqi, int size)
	{
		size_t len1 = size - 3;
		for (size_t sk = 0; sk < NAA4; sk++) taap[sk] = 0;
		for (size_t j1 = 0; j1 < len1; j1++) {
			if ((seqi[j1] >= 4) || (seqi[j1 + 1] >= 4) || (seqi[j1 + 2] >= 4) || (seqi[j1 + 3] >= 4)) continue; //skip N
			size_t c22 = seqi[j1] * NAA3 + seqi[j1 + 1] * NAA2 + seqi[j1 + 2] * NAA1 + seqi[j1 + 3];
			taap[c22]++;
		}
		for (size_t sk = 0, mm = 0; sk < NAA4; sk++) {
			aap_begin[sk] = mm;  mm += taap[sk];  taap[sk] = 0;
		}
		for (size_t j1 = 0; j1 < len1; j1++) {
			if ((seqi[j1] >= 4) || (seqi[j1 + 1] >= 4) || (seqi[j1 + 2] >= 4) || (seqi[j1 + 3] >= 4)) continue; //skip N
			size_t c22 = seqi[j1] * NAA3 + seqi[j1 + 1] * NAA2 + seqi[j1 + 2] * NAA1 + seqi[j1 + 3];
			aap_list[aap_begin[c22] + taap[c22]++] = j1;
		}
	}

};
