#!/usr/bin/env python3
# coding: utf-8

"""
Pair-wise alignment algorithm for comparison of two genome sequences

usage: python3 pairwise-alignment.py
"""

from abc import ABCMeta, abstractmethod
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        description='Pair-wise alignment algorithm for comparison of two genome sequences')
    parser.add_argument(
        '-A',
        type=int,
        default=1,
        help='score for a sequence match [1]'
    )
    parser.add_argument(
        '-B',
        type=int,
        default=1,
        help='penalty for a mismatch [1]'
    )
    parser.add_argument(
        '-D',
        type=int,
        default=2,
        help='penalty for a gap [2]'
    )
    parser.add_argument(
        '-E',
        type=int,
        default=1,
        help='affine gap penalty: a gap of size g cost -D-(g-1)*E [1]'
    )
    parser.add_argument(
        '-q', '--query',
        default='ACCGAACTACCGAAGTCGT',
        help='a query sequence ["ACCGAACTACCGAAGTCGT"]'
    )
    parser.add_argument(
        '-t', '--target',
        default='ACCGAACTGAACTCGTCGT',
        help='a target sequence ["ACCGAACTGAACTCGTCGT"]'
    )
    parser.add_argument(
        '-m', '--methods',
        choices=['nb', 'nbg', 'sw', 'swg'],
        nargs='*',
        default=['nb', 'nbg', 'sw', 'swg'],
        help='alignment methods that you want to run [nb nbg sw]'
    )
    return parser.parse_args()


class Alignment(metaclass=ABCMeta):
    def __init__(self, q, t, A, B, D):
        self.query = q
        self.target = t
        self.qlen = len(q)
        self.tlen = len(t)
        self.match_score = A
        self.mismatch_penalty = -B
        self.gap_penalty = D
        self.dp_table = self._generate_dp_table()
        self.aligned_query = ''
        self.aligned_target = ''

    @abstractmethod
    def _init_dp_table(self):
        pass

    @abstractmethod
    def _calculate_score(self):
        pass

    @abstractmethod
    def _traceback(self):
        pass

    def _generate_dp_table(self):
        return [[0 for _ in range(self.tlen + 1)] for _ in range(self.qlen + 1)]

    def _s(self, x, y):
        return self.match_score if x == y else self.mismatch_penalty

    def run(self):
        self._init_dp_table()
        self._calculate_score()
        self._traceback()
        return (self.aligned_query, self.aligned_target)


class NeedlemanWunsch(Alignment):
    def __init__(self, q, t, A, B, D):
        super().__init__(q, t, A, B, D)

    def _init_dp_table(self):
        self.dp_table[0][0] = 0
        for i in range(1, self.qlen + 1):
            self.dp_table[i][0] = -self.gap_penalty * i
        for j in range(1, self.tlen + 1):
            self.dp_table[0][j] = -self.gap_penalty * j

    def _calculate_score(self):
        for i in range(1, self.qlen + 1):
            for j in range(1, self.tlen + 1):
                self.dp_table[i][j] = max(
                    self.dp_table[i][j-1] - self.gap_penalty,
                    self.dp_table[i-1][j] - self.gap_penalty,
                    self.dp_table[i-1][j-1] + self._s(self.query[i-1], self.target[j-1])
                )

    def _traceback(self):
        i = self.qlen
        j = self.tlen
        while i > 0 and j > 0:
            if self.dp_table[i][j] == self.dp_table[i-1][j] - self.gap_penalty:
                i -= 1
                self.aligned_query += self.query[i]
                self.aligned_target += '-'
            elif self.dp_table[i][j] == self.dp_table[i][j-1] - self.gap_penalty:
                j -= 1
                self.aligned_query += '-'
                self.aligned_target += self.target[j]
            else:
                i -= 1
                j -= 1
                self.aligned_query += self.query[i]
                self.aligned_target += self.target[j]
        while i > 0:
            i -= 1
            self.aligned_query += self.query[i]
            self.aligned_target += '-'
        while j > 0:
            j -= 1
            self.aligned_query += '-'
            self.aligned_target += self.target[j]
        self.aligned_query = self.aligned_query[::-1]
        self.aligned_target = self.aligned_target[::-1]


class NeedlemanWunschGotoh(Alignment):
    def __init__(self, q, t, A, B, D, E):
        super().__init__(q, t, A, B, D)
        self.ext_penalty = E
        self.dp_table_X = self._generate_dp_table()
        self.dp_table_Y = self._generate_dp_table()

    def _init_dp_table(self):
        self.dp_table[0][0] = 0
        self.dp_table_X[0][0] = -float('inf')
        self.dp_table_Y[0][0] = -float('inf')
        for i in range(1, self.qlen + 1):
            self.dp_table[i][0] = -float('inf')
            self.dp_table_X[i][0] = -float('inf')
            self.dp_table_Y[i][0] = -self.gap_penalty - (i - 1) * self.ext_penalty
        for j in range(1, self.tlen + 1):
            self.dp_table[0][j] = -float('inf')
            self.dp_table_X[0][j] = -self.gap_penalty - (j - 1) * self.ext_penalty
            self.dp_table_Y[0][j] = -float('inf')

    def _calculate_score(self):
        for i in range(1, self.qlen + 1):
            for j in range(1, self.tlen + 1):
                self.dp_table[i][j] = max(
                    self.dp_table[i-1][j-1] + self._s(self.query[i-1], self.target[j-1]),
                    self.dp_table_X[i-1][j-1] + self._s(self.query[i-1], self.target[j-1]),
                    self.dp_table_Y[i-1][j-1] + self._s(self.query[i-1], self.target[j-1])
                )
                self.dp_table_X[i][j] = max(
                    self.dp_table[i][j-1] - self.gap_penalty,
                    self.dp_table_X[i][j-1] - self.ext_penalty
                )
                self.dp_table_Y[i][j] = max(
                    self.dp_table[i-1][j] - self.gap_penalty,
                    self.dp_table_Y[i-1][j] - self.ext_penalty
                )

    def _traceback(self):
        i = self.qlen
        j = self.tlen
        tmp_table = self._get_tmp_table()
        while i > 0 and j > 0:
            if tmp_table == 'M':
                if self.dp_table[i][j] == self.dp_table_X[i-1][j-1] + self._s(self.query[i-1], self.target[j-1]):
                    tmp_table = 'X'
                elif self.dp_table[i][j] == self.dp_table_Y[i-1][j-1] + self._s(self.query[i-1], self.target[j-1]):
                    tmp_table = 'Y'
                i -= 1
                j -= 1
                self.aligned_query += self.query[i]
                self.aligned_target += self.target[j]
            elif tmp_table == 'X':
                if self.dp_table_X[i][j] == self.dp_table[i][j-1] - self.gap_penalty:
                    tmp_table = 'M'
                j -= 1
                self.aligned_query += '-'
                self.aligned_target += self.target[j]
            else:
                if self.dp_table_Y[i][j] == self.dp_table[i-1][j] - self.gap_penalty:
                    tmp_table = 'M'
                i -= 1
                self.aligned_query += self.query[i]
                self.aligned_target += '-'
        while i > 0:
            i -= 1
            self.aligned_query += self.query[i]
            self.aligned_target += '-'
        while j > 0:
            j -= 1
            self.aligned_query += '-'
            self.aligned_target += self.target[j]
        self.aligned_query = self.aligned_query[::-1]
        self.aligned_target = self.aligned_target[::-1]

    def _get_tmp_table(self):
        max_score = max(
            self.dp_table[-1][-1], self.dp_table_X[-1][-1], self.dp_table_Y[-1][-1])
        if self.dp_table[-1][-1] == max_score:
            return 'M'
        elif self.dp_table_X[-1][-1] == max_score:
            return 'X'
        else:
            return 'Y'


class SmithWaterman(Alignment):
    def __init__(self, q, t, A, B, D):
        return super().__init__(q, t, A, B, D)

    def _init_dp_table(self):
        self.dp_table[0][0] = 0
        for i in range(1, self.qlen + 1):
            self.dp_table[i][0] = 0
        for j in range(1, self.tlen + 1):
            self.dp_table[0][j] = 0

    def _calculate_score(self):
        for i in range(1, self.qlen + 1):
            for j in range(1, self.tlen + 1):
                self.dp_table[i][j] = max(
                    0,
                    self.dp_table[i][j-1] - self.gap_penalty,
                    self.dp_table[i-1][j] - self.gap_penalty,
                    self.dp_table[i-1][j-1] + self._s(self.query[i-1], self.target[j-1])
                )
    
    def _traceback(self):
        i, j = self._get_start_subscript()
        while i > 0 and j > 0 and self.dp_table[i][j] > 0:
            if self.dp_table[i][j] == self.dp_table[i-1][j] - self.gap_penalty:
                i -= 1
                self.aligned_query += self.query[i]
                self.aligned_target += '-'
            elif self.dp_table[i][j] == self.dp_table[i][j-1] - self.gap_penalty:
                j -= 1
                self.aligned_query += '-'
                self.aligned_target += self.target[j]
            else:
                i -= 1
                j -= 1
                self.aligned_query += self.query[i]
                self.aligned_target += self.target[j]
        self.aligned_query = self.aligned_query[::-1]
        self.aligned_target = self.aligned_target[::-1]

    def _get_start_subscript(self):
        max_score = 0
        start_i = start_j = 0
        for i in range(1, self.qlen + 1):
            for j in range(1, self.tlen + 1):
                if self.dp_table[i][j] >= max_score:
                    max_score = self.dp_table[i][j]
                    start_i = i
                    start_j = j
        return start_i, start_j


class SmithWatermanGotoh(Alignment):
    def __init__(self, q, t, A, B, D, E):
        super().__init__(q, t, A, B, D)
        self.ext_penalty = E
        self.dp_table_X = self._generate_dp_table()
        self.dp_table_Y = self._generate_dp_table()

    def _init_dp_table(self):
        pass

    def _calculate_score(self):
        for i in range(1, self.qlen + 1):
            for j in range(1, self.tlen + 1):
                self.dp_table[i][j] = max(
                    0,
                    self.dp_table_X[i][j-1] + self._s(self.query[i-1], self.target[j-1]),
                    self.dp_table_Y[i-1][j] + self._s(self.query[i-1], self.target[j-1]),
                    self.dp_table[i-1][j-1] + self._s(self.query[i-1], self.target[j-1]),
                )
                self.dp_table_X[i][j] = max(
                    0,
                    self.dp_table[i][j-1] - self.gap_penalty,
                    self.dp_table_X[i][j-1] - self.ext_penalty
                )
                self.dp_table_Y[i][j] = max(
                    0,
                    self.dp_table[i-1][j] - self.gap_penalty,
                    self.dp_table_Y[i-1][j] - self.ext_penalty
                )

    def _traceback(self):
        i, j, tmp_table = self._get_start_state()
        while i > 0 and j > 0:
            if (tmp_table == 'M' and self.dp_table[i][j] <= 0) \
                    or (tmp_table == 'X' and self.dp_table_X[i][j] <= 0) \
                    or (tmp_table == 'Y' and self.dp_table_Y[i][j] <= 0):
                break
            if tmp_table == 'M':
                if self.dp_table[i][j] == self.dp_table_X[i-1][j-1] + self._s(self.query[i-1], self.target[j-1]):
                    tmp_table = 'X'
                elif self.dp_table[i][j] == self.dp_table_Y[i-1][j-1] + self._s(self.query[i-1], self.target[j-1]):
                    tmp_table = 'Y'
                i -= 1
                j -= 1
                self.aligned_query += self.query[i]
                self.aligned_target += self.target[j]
            elif tmp_table == 'X':
                if self.dp_table_X[i][j] == self.dp_table[i][j-1] - self.gap_penalty:
                    tmp_table = 'M'
                j -= 1
                self.aligned_query += '-'
                self.aligned_target += self.target[j]
            else:
                if self.dp_table_Y[i][j] == self.dp_table[i-1][j] - self.gap_penalty:
                    tmp_table = 'M'
                i -= 1
                self.aligned_query += self.query[i]
                self.aligned_target += '-'
        self.aligned_query = self.aligned_query[::-1]
        self.aligned_target = self.aligned_target[::-1]

    def _get_start_state(self):
        max_score = -1
        start_i = start_j = 0
        tmp_table = ''
        for i in range(self.qlen + 1):
            for j in range(self.tlen + 1):
                if self.dp_table[i][j] >= max_score:
                    start_i, start_j = i, j
                    max_score = self.dp_table[i][j]
                    tmp_table = 'M'
                elif self.dp_table_X[i][j] >= max_score:
                    start_i, start_j = i, j
                    max_score = self.dp_table_X[i][j]
                    tmp_table = 'X'
                elif self.dp_table_Y[i][j] >= max_score:
                    start_i, start_j = i, j
                    max_score = self.dp_table_Y[i][j]
                    tmp_table = 'Y'
        return start_i, start_j, tmp_table


def output_result(query, target):
    symbol = ''
    for q, t in zip(query, target):
        symbol += '|' if q == t and q != '-' else ' '
    target_len = len([t for t in target if t != '-'])
    match_len = len([i for i in symbol if i != ' '])
    if target_len == 0:
        print('no alignment (>_<)')
        return
    identity = match_len / target_len * 100
    print(f'query:  {query}')
    print(f'        {symbol}')
    print(f'target: {target}')
    print('')
    print(f'{identity:.2f}% ({match_len}/{target_len})')


def main():
    args = get_args()
    print('********* Input sequences **********')
    print('query:  {}'.format(args.query))
    print('target: {}'.format(args.target))
    print('')
    if 'nb' in args.methods:
        print('********* Needleman-Bunsch *********')
        q, t = NeedlemanWunsch(args.query, args.target, args.A, args.B, args.D).run()
        output_result(q, t)
        print('')
    if 'nbg' in args.methods:
        print('****** Needleman-Bunsch-Gotoh ******')
        q, t = NeedlemanWunschGotoh(args.query, args.target, args.A, args.B, args.D, args.E).run()
        output_result(q, t)
        print('')
    if 'sw' in args.methods:
        print('********** Smith-Waterman **********')
        q, t = SmithWaterman(args.query, args.target, args.A, args.B, args.D).run()
        output_result(q, t)
        print('')
    if 'swg' in args.methods:
        print('******* Smith-Waterman-Gotoh *******')
        q, t = SmithWatermanGotoh(args.query, args.target, args.A, args.B, args.D, args.E).run()
        output_result(q, t)


if __name__ == "__main__":
    main()
