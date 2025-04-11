#include <thread>
#include <vector>

void DoClustering(int T)
{
    int i, k;
    int NAA = options.NAA;
    double aa1_cutoff = options.cluster_thd;
    double aas_cutoff = 1 - (1 - options.cluster_thd) * 4;
    double aan_cutoff = 1 - (1 - options.cluster_thd) * options.NAA;
    int seq_no = sequences.size();
    int frag_no = seq_no;
    int frag_size = options.frag_size;
    std::valarray<size_t>  letters(T);

    if (frag_size) {
        frag_no = 0;
        for (i = 0; i < seq_no; i++) frag_no += (sequences[i]->size - NAA) / frag_size + 1;
    }

    if (not options.isEST)
        cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd, options.tolerance, NAA);

    std::vector<WorkingParam> params(T);
    std::vector<WorkingBuffer> buffers(T);
    for (i = 0; i < T; i++) {
        params[i] = WorkingParam(aa1_cutoff, aas_cutoff, aan_cutoff);
        buffers[i].Set(frag_no, max_len);
    }

    WordTable word_table(options.NAA, NAAN);
    WordTable last_table(options.NAA, NAAN);

    int N = sequences.size();
    size_t mem_need = MinimalMemory(frag_no, buffers[0].total_bytes, T);
    size_t mem;
    size_t tabsize = 0;
    int remaining = 0;

    options.ComputeTableLimits(min_len, max_len, len_n50, mem_need);

    for (i = 0; i < N; ) {
        int start = i;
        int m = i;
        size_t sum = remaining;
        float redundancy = (rep_seqs.size() + 1.0) / (i + 1.0);
        size_t max_items = options.max_entries;
        size_t max_seqs = options.max_sequences;
        size_t items = 0;
        if (i == 0 && max_seqs > 1000) {
            max_items /= 8;
            max_seqs /= 8;
        }
        while (m < N && (sum * redundancy) < max_seqs && items < max_items) {
            Sequence* seq = sequences[m];
            if (!(seq->state & IS_REDUNDANT)) {
                if (options.store_disk) seq->SwapIn();
                items += static_cast<size_t>(seq->size * redundancy);
                sum += 1;
            }
            m++;
        }
        if ((m > i + 1E4) && (m > i + (N - i) / (2 + T))) m = i + (N - i) / (2 + T);
        if (m == i || m >= N) {
            m = N;
            if (m > i + 1E3) m = i + (N - i) / (2 + T);
        }
        printf("\r# comparing sequences from  %9i  to  %9i\n", i, m);
        if (last_table.size) {
            int print = (m - i) / 20 + 1;
            std::vector<std::thread> threads;
            for (int j = i; j < m; j++) {
                threads.emplace_back([&, j]() {
                    Sequence* seq = sequences[j];
                    if (seq->state & IS_REDUNDANT) return;
                    int tid = std::this_thread::get_id().hash() % T;
                    CheckOne(seq, last_table, params[tid], buffers[tid]);
                    if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
                    if (j % print == 0) {
                        std::cout << ".";
                        std::cout.flush();
                    }
                });
            }
            for (auto& thread : threads) {
                thread.join();
            }
            int may_stop = 0;
            int self_stop = 0;
            float p0 = 0;
            int min = last_table.sequences[last_table.sequences.size() - 1]->size;
            int m0 = m;
            bool stop = false;
            std::vector<std::thread> threads2;
            for (int j = m - 1; j < N; j++) {
                threads2.emplace_back([&, j]() {
                    if (!stop) {
                        if (j + 1 == N) may_stop = 1;
                        if (j == (m0 - 1)) {
                            int tid = std::this_thread::get_id().hash() % T;
                            for (int ks = i; ks < m; ks++) {
                                Sequence* seq = sequences[ks];
                                i = ks + 1;
                                if (seq->state & IS_REDUNDANT) continue;
                                ClusterOne(seq, ks, word_table, params[tid], buffers[tid]);
                                if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
                                if (may_stop and word_table.sequences.size() >= 100) break;
                                if (word_table.size >= max_items) break;
                                int tmax = max_seqs - (frag_size ? seq->size / frag_size + 1 : 0);
                                if (word_table.sequences.size() >= tmax) break;
                            }
                            self_stop = 1;
                        } else {
                            Sequence* seq = sequences[j];
                            if (seq->state & IS_REDUNDANT) continue;
                            if (options.store_disk) {
                                seq->SwapIn();
                            }
                            int tid = std::this_thread::get_id().hash() % T;
                            CheckOne(seq, last_table, params[tid], buffers[tid]);
                            if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
                            if (min > params[tid].len_upper_bound) {
                                may_stop = 1;
                                stop = true;
                            }
                            if (self_stop && tid == 1) {
                                float p = (100.0 * j) / N;
                                if (p > p0 + 1E-1) {
                                    printf("\r%4.1f%%", p);
                                    p0 = p;
                                }
                            }
                        }
                    }
                });
            }
            for (auto& thread : threads2) {
                thread.join();
            }
        }
        if (i == start || m == N) {
            for (k = i; k < m; ) {
                int mm = k;
                auto sum = 0;
                while (mm < m && sum < 1E5) {
                    if (!(sequences[mm]->state & IS_REDUNDANT)) sum += sequences[mm]->size;
                    mm += 1;
                }
                if (mm < k + 1000) mm = k + 1000;
                if (mm > m) mm = m;
                std::vector<std::thread> threads3;
                for (auto kk = k; kk < mm; kk++) {
                    threads3.emplace_back([&, kk]() {
                        Sequence* seq = sequences[kk];
                        if (seq->state & IS_REDUNDANT) return;
                        int tid = std::this_thread::get_id().hash() % T;
                        CheckOne(seq, word_table, params[tid], buffers[tid]);
                        if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
                    });
                }
                for (auto& thread : threads3) {
                    thread.join();
                }
                bool bk = false;
                for (int ks = k; ks < mm; ks++) {
                    Sequence* seq = sequences[ks];
                    i = k = ks + 1;
                    if (seq->state & IS_REDUNDANT) continue;
                    ClusterOne(seq, ks, word_table, params[0], buffers[0]);
                    bk = true;
                    if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
                    if (word_table.size >= max_items) break;
                    int tmax = max_seqs - (frag_size ? seq->size / frag_size + 1 : 0);
                    if (word_table.sequences.size() >= tmax) break;
                    bk = false;
                }
                if (bk) break;
            }
        } else if (i < m) {
            remaining = remaining / 2 + (m - i);
            printf("\r---------- %6i remaining sequences to the next cycle\n", m - i);
        }
        std::println("---------- new table with {:8} representatives", word_table.sequences.size());
        if ((last_table.size + word_table.size) > tabsize)
            tabsize = last_table.size + word_table.size;
        last_table.Clear();
        last_table.sequences.swap(word_table.sequences);
        last_table.indexCounts.swap(word_table.indexCounts);
        last_table.size = word_table.size;
        word_table.size = 0;
    }
    std::println("\n{:9}  finished  {:9}  clusters", sequences.size(), rep_seqs.size());
    mem = (mem_need + tabsize * sizeof(IndexCount)) / MEGA_MiBi;
    std::println("\nApproximated maximum memory consumption: {}M", mem);
    last_table.Clear();
    word_table.Clear();
}
