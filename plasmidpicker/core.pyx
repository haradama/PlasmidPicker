# coding: utf-8

from mmh3 import hash128
import numpy as np
cimport numpy as np
import itertools

DTYPE_INT = np.int
DTYPE_FLOAT = np.float

ctypedef np.int_t DTYPE_INT_t
ctypedef np.float_t DTYPE_FLOAT_t

cdef class CGR:
    """
    
    """
    cdef int k, array_size, num_kmers
    cdef np.ndarray cgr_array
    def __cinit__(self, int k_length):        
        self.k = k_length
        self.array_size = int((4**self.k)**0.5)
        self.cgr_array = np.zeros((self.array_size, self.array_size), dtype=DTYPE_INT)
        self.num_kmers = 0
    
    def __len__(self):
        return self.num_kmers
    
    cpdef void add_kmer(self, str kmer):
        assert len(kmer) == self.k, "The added kmer is not a {0}-mer but a {1}-mer".format(str(self.k), str(len(kmer)))
        
        cdef str nucl
        cdef int index, maxX, maxY, locX, locY
        if 'N' not in kmer:
            maxX, maxY = self.array_size, self.array_size
            locX, locY = 1, 1

            for nucl in kmer:
                if nucl == "T":
                    locX += maxX // 2
                elif nucl == "C":
                    locY += maxY // 2
                elif nucl == "G":
                    locX += maxX // 2
                    locY += maxY // 2

                maxX //= 2
                maxY //= 2

            self.cgr_array[locY-1, locX-1] += 1
            maxX, maxY = self.array_size, self.array_size
            locX, locY = 1, 1
            self.num_kmers += 1
            
    cpdef np.ndarray[DTYPE_FLOAT_t, ndim=2] get_cgr_array(self):
        if self.num_kmers == 0:
            return self.cgr_array
        else:
            return self.cgr_array / self.num_kmers
    
    cpdef void add_seq(self, str seq):
        cdef int index
        self.cgr_array = np.zeros((self.array_size, self.array_size), dtype=DTYPE_INT)
        self.num_kmers = 0
        for index in range(len(seq) - (self.k - 1)):
            self.add_kmer(seq[index:index + self.k])

    cpdef get_signatures_table(self):
        cdef int maxX, maxY, locX, locY
        cdef str nucl
        signatures_table = np.full((self.array_size, self.array_size), "N" * self.k)
        base = ("A", "T", "G", "C")
        maxX, maxY = self.array_size, self.array_size
        locX, locY = 1, 1

        for p in itertools.product(*([base] * self.k)):
            for nucl in p:
                if nucl == "T":
                    locX += maxX // 2
                elif nucl == "C":
                    locY += maxY // 2
                elif nucl == "G":
                    locX += maxX // 2
                    locY += maxY // 2

                maxX //= 2
                maxY //= 2

            signatures_table[locY-1, locX-1] = "".join(p)
            maxX, maxY = self.array_size, self.array_size
            locX, locY = 1, 1

        return signatures_table

cdef class DenoisedCGR(CGR):
    cpdef set kmer_set
    def __cinit__(self, int k_length):
        CGR.__init__(k_length)
        self.kmer_set = set()
        
    cpdef void add(self, str kmer):
        assert len(kmer) == self.k, "The added kmer is not a {0}-mer but a {1}-mer".format(str(self.k), str(len(kmer)))
        if kmer in self.kmer_set:
            self.add_kmer(kmer)
        else:
            self.kmer_set.add(kmer)
            
    cpdef void add_seq(self, str seq):
        cdef int index
        self.cgr_array = np.zeros((self.array_size, self.array_size), dtype=DTYPE_INT)
        for index in range(len(seq) - (self.k - 1)):
            self.add(seq[index:index + self.k])

def kmerIter(str seq, int k_length):
    cdef str fragment
    cdef int i
    for i in range(len(seq) - (k_length - 1)):
        fragment = seq[i:i+k_length]
        yield fragment

cdef class MinHash:
    cdef int num, k_length
    def __cinit__(self, int _num, int _k_length):
        self.num = _num
        self.k_length = _k_length
    
    cpdef np.ndarray[double, ndim=2] getMin(self, str seq):
        cdef str kmer
        cdef int i, j, kmer_num
        kmer_num = len(seq) - (self.k_length - 1)
        cdef np.ndarray[double, ndim=2] hash_matrix = np.empty((self.num, kmer_num), dtype=np.double)
        
        for i in range(kmer_num):
            kmer = seq[i:i+self.k_length]
            for j in range(self.num):
                hash_matrix[j, i] = hash128(kmer, j)
        
        return hash_matrix.min(axis=1)
    
    cpdef double get_concordance(self, np.ndarray[double, ndim=1] array1, np.ndarray[double, ndim=1] array2):
        return (array1 == array2).sum() / self.num

cpdef double get_GC_content(str seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

cpdef str get_reverse_complement(str seq):
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
