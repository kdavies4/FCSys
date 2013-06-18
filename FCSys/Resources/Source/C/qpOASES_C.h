/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2009 by Hans Joachim Ferreau et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file INTERFACES/C/qpOASES_C.h
 *	\author Joel Andersson, Hans Joachim Ferreau
 *	\version 2.x
 *	\date 2009
 *
 *
 */

#ifndef QPOASES_QPOASES_C_H
#define QPOASES_QPOASES_C_H

void* qpoases_create(   char *mode, 
			int nV,
			int nC);

void* qpoases_createMPC(
			int y_hor, // Prediction horizon
			int u_hor, // Control horizon
			int nx,    // Number of states
			int nu,    // Number of inputs
			int ny,    // Number of outputs
			
			double *A,      double *B,      double *C,
			double *A_dist, double *B_dist, double *C_dist,
			
 			double *op_states,
			double *op_inputs,
			double *op_outputs,
			double *setpoints,

			double *u_min, double *u_max);

int qpOASES_feedbackMPC(
			void *object,
			double *x,
			double *u_old,
			double *delta_u);

int qpoases_init(	void *object, 
			double* H, 
			double* g, double* A, 
			double* lb, double* ub, double* lbA, double* ubA,
			int nV, int nC, int nWSRmax);

int qpoases_hotstart(	void *object,
			double* H, 
			double* g, double* A, 
			double* lb, double* ub, double* lbA, double* ubA);

int qpoases_getInfo(    void *object,
			double* obj, double* x, double* y, int* nWSR);

int qpoases_free(       void *object);

int qpoases_initMPC(    void *object,
			double *A, double *B,
			double *C, double *D,
			int np, int nu,
			double *op_states,
			double *op_inputs,
			double *op_outputs,
			double *setpoints);

#endif // QPOASES_QPOASES_C_H