export module WordTable;

import std;
import common;
import Sequence;
import Options;
import ScoreMatrix;

export struct WorkingBuffer
{
	VectorInt taap;
	VectorInt word_encodes;
	VectorIntX word_encodes_backup;
	VectorIntX word_encodes_no;
	VectorIntX aap_list;
	VectorIntX aap_begin;
	VectorInt32 indexMapping;
	MatrixInt64 score_mat;
	MatrixInt back_mat;
	VectorInt diag_score;
	VectorInt diag_score2;
	VectorInt aan_list_comp;
	std::string seqi_comp{};
	std::vector<IndexCount> lookCounts;
	size_t total_bytes;

	explicit WorkingBuffer(size_t frag = 0, size_t maxlen = 0)
	{
		Set(frag, maxlen);
		seqi_comp.resize(MAX_SEQ);
	}
	void Set(size_t _frag, size_t max_len)
	{
	    constexpr auto MAX_UAA_SQ = MAX_UAA * MAX_UAA;
	    constexpr auto MAX_UAA_SQSQ = MAX_UAA_SQ * MAX_UAA_SQ;
		const auto UAA_SQUARED = options.isEST ? MAX_UAA_SQSQ : MAX_UAA_SQ;
        const auto band = std::min<size_t>(max_len * max_len, options.band_width);
		/* each table can not contain more than MAX_TABLE_SEQ representatives or fragments! */
	    const auto frag = std::min<size_t>(_frag, MAX_TABLE_SEQ);
        taap.resize(UAA_SQUARED);
		aap_list.resize(max_len);
		aap_begin.resize(UAA_SQUARED);
		word_encodes.resize(max_len);
		word_encodes_no.resize(max_len);
		word_encodes_backup.resize(max_len);
		lookCounts.resize(frag + 2);
		indexMapping.resize(frag + 2);
		diag_score.resize(StaticOptions::MAX_DIAG);
		diag_score2.resize(StaticOptions::MAX_DIAG);
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
		total_bytes += indexMapping.size() * sizeof(std::uint32_t);
		// total_bytes += indexCounts.size()*sizeof(IndexCount);
		total_bytes += lookCounts.size() * sizeof(IndexCount);
		total_bytes += max_len * (band * sizeof(int) + sizeof(VectorInt));
		total_bytes += max_len * (band * sizeof(int) + sizeof(VectorInt64));
	}

	int EncodeWords(const Sequence& seq, int NAA, bool est)
	{
		const auto& sseq = seq.get_data();
		const char* seqi = sseq.data();
		size_t len = sseq.size();
		// check_word_encodes
		auto aan_no = len - NAA + 1;
		auto skip = 0;
		for (size_t j = 0; j < aan_no; j++) {
			const char* word = seqi + j;
			int encode = 0;
			for (auto k = 0, k1 = NAA - 1; k < NAA; k++, k1--) encode += word[k] * NAAN_array[k1];
			word_encodes[j] = word_encodes_backup[j] = encode;
		}

		if (est) {
			for (auto j = 0; j < len; j++) {
				if (seqi[j] >= 4) {                      // here N is 4
					auto i0 = (j - NAA + 1 > 0) ? j - NAA + 1 : 0;
					auto i1 = j < aan_no ? j : aan_no - 1;
					for (auto i = i0; i <= i1; i++) word_encodes[i] = -1;
				}
			}
			for (auto j = 0; j < aan_no; j++) skip += (word_encodes[j] == -1);
		}

		std::sort(word_encodes.begin(), word_encodes.begin() + aan_no);
		for (auto j = 0; j < aan_no; j++) word_encodes_no[j] = 1;
		for (auto j = aan_no - 1; j; j--) {
			if (word_encodes[j] == word_encodes[j - 1]) {
				word_encodes_no[j - 1] += word_encodes_no[j];
				word_encodes_no[j] = 0;
			}
		}
		return skip;
		// END check_word_encodes
	}
	void ComputeAAP(const std::string_view seqi)
	{
		size_t len1 = seqi.size();
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
	void ComputeAAP2(const std::string_view seqi)
	{
		size_t len1 = seqi.size() - 2;
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


export struct WordTable {
    const size_t NAA;  // length of word
    const size_t NAAN; // rows of table
    const size_t frag_size;

    std::vector<NVectorIndexCount> indexCounts; // hold index and word counts of seqs
    std::vector<Sequence*> sequences;
    size_t size = 0;
    size_t frag_count = 0;


	WordTable(size_t naa, size_t naan): frag_size(options.frag_size), NAA(naa), NAAN(naan), indexCounts(naan) {}

    void Clear() {
        size = 0;
        frag_count = 0;
        sequences.clear();
        indexCounts.clear();
    }

    size_t AddWordCounts(const NVectorIndexCount& counts, Sequence& seq, bool skipN = false) {
        auto aan_no = counts.size();
        auto idx = sequences.size();
        for (auto i = 0; i < aan_no; i++) {
            IndexCount ic = counts[i];
            auto k = ic.count;
            if (k) {
                auto j = ic.index;
                if (skipN && j < 0)
                    continue; // for those has 'N'
                NVectorIndexCount row = indexCounts[j];
                ic.index = idx;
                row.push_back(ic);
                size += 1;
            }
        }
        sequences.push_back(&seq);
        return OK_FUNC;
    }

    size_t AddWordCounts(const size_t aan_no, const VectorInt& word_encodes, VectorIntX& word_encodes_no, size_t idx) {
        int i, j, k;
        const bool skipN = options.isEST;
        for (i = 0; i < aan_no; i++) {
            if ((k = word_encodes_no[i])) {
                j = word_encodes[i];
                if (skipN && j < 0)
                    continue; // for those has 'N'
                NVectorIndexCount& row = indexCounts[j];
                row.push_back(IndexCount(idx, k));
                size += 1;
                // if( k >1 ) printf( " %3i", k );
            }
        }
        // printf( "\n" );
        return OK_FUNC;
    }

    size_t AddWordCountsFrag(const size_t aan_no, VectorIntX& word_encodes, VectorIntX& word_encodes_no, size_t num_frags) {
	    auto word_encode_iter = VectorIterFactory{word_encodes, frag_size};
        for (size_t i = 0; i < num_frags; ++i) {
            auto end_i = i + 1;
            auto frag_end = end_i * frag_size < word_encodes.size() ? end_i * (frag_size - 1) : (aan_no - 1);
            std::sort(word_encode_iter(i * frag_size), word_encode_iter(frag_end + 1));
        }

        for (size_t j = aan_no - 1; j; j--) {
            if (word_encodes[j] == word_encodes[j - 1]) {
                word_encodes_no[j - 1] += word_encodes_no[j];
                word_encodes_no[j] = 0;
            }
        }
        // END check_word_encodes

        for (size_t i = 0; i < aan_no; i += frag_size) {
            auto k = frag_size < (aan_no - i) ? frag_size : (aan_no - i);
            auto fra = i / frag_size;
            // AddWordCounts(k, word_encodes+i, word_encodes_no+i, NR90f_no+fra);
            for (auto i1 = i; i1 < i + k; i1++) {
                if (auto k1 = word_encodes_no[i1]) {
                    auto j = word_encodes[i1];
                    NVectorIndexCount& row = indexCounts[j];
                    row.emplace_back(IndexCount(frag_count + fra, k1));
                    size += 1;
                }
            }
        }
        frag_count += num_frags;

        return 0;
    }

int CountWords(int total_words, WorkingBuffer& buf, bool is_est = false, int min_threshold = 0) {
	    // Reset index mapping for all indices in lookup_counts
	    for (const auto& index_count : buf.lookCounts)
	        buf.indexMapping[index_count.index] = 0;
	    buf.lookCounts.clear();

	    // Start processing from the first word
	    auto& current_word = buf.word_encodes[0];
	    auto current_index = 0;

	    // Skip words with 'N' if is_est is true
	    if (is_est)
	        while (current_word < 0)
	            current_index++, current_word++;

	    auto& current_word_count = buf.word_encodes[current_index];

	    // Iterate through all words
	    for (; current_index < total_words; current_index++, current_word++, current_word_count++) {
	        auto word_index = current_word;
	        auto word_count = current_word_count;

	        // Skip words with zero count
	        if (word_count == 0)
	            continue;

	        // Retrieve the list of index-count pairs for the current word
	        NVectorIndexCount& word_index_counts = indexCounts[word_index];
	        auto word_index_count_size = word_index_counts.size();
	        int remaining_words = total_words - current_index + 1;

	        // Process each index-count pair for the current word
	        for (const auto& index_count_pair : word_index_counts) {
	            size_t min_count = index_count_pair.count < word_count ? index_count_pair.count : word_count;

	            //     IndexCount *ic = buf.lookCounts.data();
	            //     for(j=0; j<buf.lookCounts.size; j++, ic++) buf.indexMapping[ ic->index ] = 0;
	            //     buf.lookCounts.size = 0;

	            // Resize index mapping if necessary
	            if (buf.indexMapping.size() <= index_count_pair.index)
	                buf.indexMapping.resize(index_count_pair.index + 1, 0);

	            auto& mapping_value = buf.indexMapping[index_count_pair.index];

	            // If the index is not yet mapped
	            if (mapping_value == 0) {
	                // Skip if the remaining words are below the threshold
	                if (remaining_words < min_threshold)
	                    continue;

	                // Add a new entry to lookup_counts
	                buf.lookCounts.emplace_back(index_count_pair.index, min_count);
	                mapping_value = buf.lookCounts.size();
	            } else {
	                // Update the count for an existing entry
	                // indexes are 1-based in lookCounts
	                buf.lookCounts[mapping_value - 1].count += min_count;
	            }
	        }
	    }

	    // Mark the end of lookup_counts
	    buf.lookCounts[buf.lookCounts.size()].count = 0;

	    return OK_FUNC;
	}
    void PrintAll() {
        int i, j, k;
        int cols = 0;
        long long total_words = 0;
        k = 0;
        for (i = 0; i < NAAN; i++) {
            int size = indexCounts[i].size();
            if (size == 0)
                continue;
            cols++;
            std::cout << k << "\t" << i << "\tsize:" << size << "\t";
            for (j = 0; j < size; j++) {
                std::cout << indexCounts[i][j].index << "," << indexCounts[i][j].count << " ";
                total_words += indexCounts[i][j].count;
            }
            std::cout << std::endl;
            k++;
        }

        std::cout << "total cols: " << cols << " total words: " << total_words << std::endl;
    }
};

// WordTable::CountWords(int aan_no, WorkingBuffer& buf, bool est, int min)
// {
//     // int S = frag_count ? frag_count : sequences.size();
//     int  j, k, j0, j1, k1, m;
//     int ix1, ix2, ix3, ix4;
//     IndexCount tmp;
//
//     IndexCount *ic = buf.lookCounts.data();
//     for(j=0; j<buf.lookCounts.size; j++, ic++) buf.indexMapping[ ic->index ] = 0;
//     buf.lookCounts.size = 0;
//
//     int *we = & buf.word_encodes[0];
//     j0 = 0;
//     if( est ) while( *we <0 ) j0++, we++; // if met short word has 'N'
//     INTs *wen = & buf.word_encodes_no[j0];
//     //printf( "\nquery : " );
//     for (; j0<aan_no; j0++, we++, wen++) {
//         j  = *we;
//         j1 = *wen;
//         //if( j1 >1 ) printf( " %3i", j1 );
//         if( j1==0 ) continue;
//         auto& one = indexCounts[j];
//         k1 = one.Size();
//         IndexCount *ic = one.items;
//
//         int rest = aan_no - j0 + 1;
//         for (k=0; k<k1; k++, ic++){
//             int c = ic->count < j1 ? ic->count : j1;
//             size_t *idm = buf.indexMapping.data() + ic->index;
//             if( *idm ==0 ){
//                 if( rest < min ) continue;
//                 IndexCount *ic2 = buf.lookCounts.items + buf.lookCounts.size;
//                 buf.lookCounts.size += 1;
//                 *idm = buf.lookCounts.size;
//                 ic2->index = ic->index;
//                 ic2->count = c;
//             }else{
//                 buf.lookCounts[ *idm - 1 ].count += c;
//             }
//         }
//     }
//     //printf( "%6i %6i\n", S, lookCounts.size );
//     buf.lookCounts[ buf.lookCounts.size ].count = 0;
//     //printf( "\n\n" );
//     return OK_FUNC;
// }
