/*
 * Copyright 2014 Christopher Nilsen
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef LAPACK_H
#define LAPACK_H

void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr,
		double *wi, double *vl, int *ldvl, double *vr, int *ldvr,
		double *work, int *lwork, int *info);
void dlasrt_(char *id, int *n, double *x, int *info);
void dsygv_(int *itype, char *jobz, char *uplo, int *n, double *a, int *lda,
		double *b, int *ldb, double *w, double *work, int *lwork,
		int *info);

#endif
