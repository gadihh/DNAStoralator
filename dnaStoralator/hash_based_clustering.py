import argparse
import os
import shutil
import threading

from simulator import *
import time
import math

class HashBasedCluster:

    def __init__(self, chosen_technology):
        self.technology = chosen_technology
        if platform.system() == "Linux":
            self.shuffled_file = '/home_nfs/sgamer8/DNAindex' + str(
                self.index) + '/files/' + self.technology + '/' + 'errors_shuffled.txt'
        elif platform.system() == "Windows":
            self.shuffled_file = 'files/' + self.technology + '/' + 'errors_shuffled.txt'

        if platform.system() == "Linux":
            self.evyat_path = '/home_nfs/sgamer8/DNAindex' + str(self.index) + '/files/' + self.technology + '/' + 'evyat.txt'
        elif platform.system() == "Windows":
            self.evyat_path = 'files/' + self.technology + '/' + 'evyat.txt'

    def hash_fun(self, x, a, w, l):
        ind = x.find(a)
        return x[ind:min(len(x), ind + w + l)]

    def ind_st(self, st):
        # All 3-grams: AAA,AAG,AAC,AAT,AGA,AGC,...,TTA,TTG,TTC,TTT
        # A-0, G-1, C-2, T-3
        # index of CAT = Decimal representaion of (203)
        N_q = {"A": 0, "C": 1, "G": 2, "T": 3}
        dec = 0
        for i in range(0, len(st)):
            dec += N_q[st[i]] * (4 ** (len(st) - i - 1))
        return dec

    def bin_sig(self, x, q):
        bs = [0] * (4 ** q)
        for i in range(0, len(x) - q + 1):
            st = x[i:i + q]
            bs[self.ind_st(st)] = 1
        bs_str = ''.join(str(e) for e in bs)
        return bs_str

    def ham_dis(self, x, y):
        dis = 0
        for i in range(0, len(x)):
            if x[i] != y[i]:
                dis += 1
        return dis

    def rand_perm(self, w):
        return ''.join(random.choice('ACGT') for _ in range(w))

    def rep_find(self, inp, parent):
        temp = inp
        cnt = 0
        while parent[temp] != temp and cnt < 10:
            cnt += 1
            temp = parent[temp]
        return temp

    def edit_dis(self, s1, s2):
        m = len(s1) + 1
        n = len(s2) + 1

        tbl = {}
        for i in range(m):
            tbl[i, 0] = i
        for j in range(n):
            tbl[0, j] = j
        for i in range(1, m):
            for j in range(1, n):
                cost = 0 if s1[i - 1] == s2[j - 1] else 1
                tbl[i, j] = min(tbl[i, j - 1] + 1, tbl[i - 1, j] + 1, tbl[i - 1, j - 1] + cost)

        return tbl[i, j]

    def display_parent(self, parent):
        clstr = {}
        for i in range(0, len(parent)):
            clstr[i] = []
        for i in range(0, len(parent)):
            clstr[self.rep_find(i, parent)].append(i)
        return clstr

    def min_max(self, val1, val2):
        min_val = min(val1, val2)
        max_val = max(val1, val2)
        return min_val, max_val

    def hash_based_cluster(self, number_of_steps, windows_size, similarity_threshold, report_func):
        reads_err = []
        strand_id = 0
        original_strand_dict = {} #map from orig strand id to the actual strand
        reads_err_original_strand_dict = {}  # map from read_err to it's orig strand id
        temp_evyat_path = self.evyat_path.removesuffix('evyat.txt')
        temp_evyat_path += 'temp_evyat.txt'
        with open(self.evyat_path, 'r') as evyat_f:
            line = "j"
            while line:
                line = evyat_f.readline()
                if not line:
                    break
                original_strand_dict.update({strand_id: line})
                line = evyat_f.readline()  # line == '*****************************\n':
                line = evyat_f.readline()
                while line != '\n':
                    reads_err.append(line.strip())
                    reads_err_original_strand_dict.update({len(reads_err) : strand_id})
                    line = evyat_f.readline()
                strand_id = strand_id+1
                line = evyat_f.readline()


        read_err_dict = {}
        reads_err_ind = [0] * (len(reads_err))
        parent = [0] * (len(reads_err))
        bin_sig_arr = []

        local_comm = number_of_steps
        tmp_bar = int(local_comm / 10)
        total_bar_size = local_comm + 2 * tmp_bar

        q = windows_size
        for i in range(0, len(reads_err)):
            if i % 100 == 0:
                report_func(total_bar_size, 1 + int((tmp_bar * i) / len(reads_err)))
            read_err_dict.update({i: reads_err[i]})
            reads_err_ind[i] = (i, reads_err[i])
            parent[i] = i
            bin_sig_arr.append(self.bin_sig(reads_err[i], q))

        C_til = self.display_parent(parent)

        dist_arr = []
        # th_low, th_high:
        for i in range(1, len(reads_err)):
            if i % 100 == 0:
                report_func(total_bar_size, tmp_bar + int((tmp_bar * i) / len(reads_err)))
            dist_arr.append(self.ham_dis(bin_sig_arr[i], bin_sig_arr[0]))
        dist_arr.sort()
        for i in range(0, 1000):
            if dist_arr[i + 1] - dist_arr[i] > 10:
                break
        th_low = min(dist_arr[i] + 5, dist_arr[i + 1] - 8)
        th_high = min(dist_arr[i + 1] - 5, th_low + 10)

        r = similarity_threshold
        w = math.ceil(math.log(len(reads_err[0]), 4))
        l = math.ceil(math.log(len(reads_err), 4))
        for lcl_step in range(0, local_comm):
            report_func(total_bar_size, (2 * tmp_bar) + lcl_step)
            hash_select = [0] * (len(reads_err))
            for i in range(0, len(C_til)):
                if len(C_til[i]) != 0:
                    hash_select[random.choice(C_til[i])] = 1

            a = self.rand_perm(w)
            hash_C_til = [0] * (len(reads_err))
            for i in range(0, len(reads_err_ind)):
                if hash_select[i] == 0:
                    hash_C_til[i] = (i, "")
                else:
                    hash_C_til[i] = (i, self.hash_fun(reads_err_ind[i][1], a, w, l))
            hash_C_til.sort(key=lambda x: x[1])

            cnt = 0
            for i in range(0, len(hash_C_til) - 1):
                if hash_C_til[i][1] == "":
                    continue
                else:
                    if hash_C_til[i][1] == hash_C_til[i + 1][1]:
                        x = reads_err[hash_C_til[i][0]]
                        y = reads_err[hash_C_til[i + 1][0]]
                        if ((self.ham_dis(bin_sig_arr[hash_C_til[i][0]], bin_sig_arr[hash_C_til[i + 1][0]]) <= th_low) or
                                ((self.ham_dis(bin_sig_arr[hash_C_til[i][0]],
                                               bin_sig_arr[hash_C_til[i + 1][0]]) <= th_high) and
                                 self.edit_dis(x, y) <= r)):
                            cnt += 1
                            min_temp, max_temp = self.min_max(self.rep_find(hash_C_til[i][0], parent),
                                                              self.rep_find(hash_C_til[i + 1][0], parent))
                            C_til[min_temp].extend(C_til[max_temp])
                            C_til[max_temp] = []
                            parent[max_temp] = min_temp

        clusters = [sorted(x) for x in list(C_til.values()) if x != []]
        with open(temp_evyat_path, 'w', newline='\n') as temp_f:
            for cluster in clusters:
                orig_strand_candidates = []
                for cluster_element in cluster:
                    orig_strand_candidates.append(reads_err_original_strand_dict.get(cluster_element))
                orig_strand_id = max(orig_strand_candidates, key= orig_strand_candidates.count)
                temp_f.write(str(original_strand_dict.get(orig_strand_id)))
                temp_f.write('*****************************\n')
                for cluster_element in cluster:
                    temp_f.write(read_err_dict.get(cluster_element) + '\n')
                temp_f.write('\n\n')

        os.remove(self.evyat_path)
        os.rename(temp_evyat_path, self.evyat_path)