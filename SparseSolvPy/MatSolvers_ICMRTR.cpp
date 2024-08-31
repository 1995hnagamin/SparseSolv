
#include "MatSolvers.hpp"
#include "SparseMatOperators.hpp"

/* ��p���O��� */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
MRTR�\���o
//=======================================================
//=======================================================
*/


/*//=======================================================
// �� MRTR�ŉ���(�E�Ӄm���������v�Z�p�^�[��)
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const double *vecB, double *results, bool init){
	const slv_int size = size0;
	double norm = 0;
	for(int i = 0 ; i < size ; i++){
		norm += vecB[i]*vecB[i];
	}
	norm = sqrt(norm);
	return solveICMRTR(size, conv_cri, max_ite, accera, norm, matA, vecB, results, init);
}


/*//=======================================================
// �� MRTR�ŉ���(���͉E�ӂ�Vector)
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const std::vector<double>& vecB, std::vector<double>& results, bool init){
	double* vecBa = new double[size0];
	double* results_a = new double[size0];
	double norm=0;
	for(int i = 0 ; i < size0 ; i++){
		const double temp = vecB[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = results[i];
	}
	norm = sqrt(norm);
	bool bl = solveICMRTR(size0, conv_cri, max_ite, accera, norm, matA, vecBa, results_a, init);
	for(int i = 0 ; i < size0 ; i++){
		results[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}

/*//=======================================================
// �� MRTR�ŉ���(���͉E�ӂ�Eigen)
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init){
	double* vecBa = new double[size0];
	double norm=0;
	for(int i = 0 ; i < size0 ; i++){
		const double temp = vecB(i);
		vecBa[i] = temp;
		norm += temp*temp;
	}
	norm = sqrt(norm);
	bool bl = solveICMRTR(size0, conv_cri, max_ite, accera, norm, matA, vecBa, results, init);
	delete[] vecBa;
	return bl;
}

/*//=======================================================
// �� MRTR�ŉ����E�O�����s�{��
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init){
	/* �Ίp�X�P�[�����O���� */
	if(is_diag_scale){
		return solveICMRTR_diag(size0, conv_cri, max_ite, accera, normB, matA, vecB, results, init);
	}

	/* �R���X�L�[�p�X�p�[�X�s��쐬 */
	double *diagD = new double[size0];
	double accela_val = accera;
	SparseMat matL;
	/* �����W�������Ȃ玩�����胂�[�h�� */
	if(accera < -1){
		accela_val = 1.05;
		/* �Ίp�����ɂȂ�܂Ŏ��{ */
		for(int kkk = 0; kkk < 10; kkk++){
			matL = matA.IC_decomp(diagD, accela_val);
			bool ok = true;
			for(slv_int i = 0; i < size0; i++){
				ok &= (diagD[i] > 0);
			}
			if(ok){
				break;
			}
			accela_val += 0.05;
		}
	}else{
		matL = matA.IC_decomp(diagD, accera);
	}
	/* �m�� */
	SparseMat matL_tr = matL.trans();

	bool bl= solveICMRTR(size0, conv_cri, max_ite, normB, diagD, matA.matrix, matL.matrix, matL_tr.matrix, vecB, results, init);
	delete[] diagD;
	return bl;	
}

/*//=======================================================
// �� MRTR�ŉ����E�O�����s�{��
//=======================================================*/
bool MatSolvers::solveICMRTR_diag(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init){
	/* �Ίp�X�P�[�����O */
	double* vecB2 = new double[size0];
	SparseMat matD = matA.diagScaling(vecB2, vecB);
	SparseMat matDAD = matD*matA*matD;

	double normB2 = 0;
	for(int i = 0; i < size0; i++){
		normB2 += vecB2[i]*vecB2[i];
	}
	normB2 = sqrt(normB2);

	/* �R���X�L�[�p�X�p�[�X�s��쐬 */
	double* diagD = new double[size0];
	double accela_val = accera;
	SparseMat matL;
	/* �����W�������Ȃ玩�����胂�[�h�� */
	if(accera < -1){
		accela_val = 1.05;
		/* �Ίp�����ɂȂ�܂Ŏ��{ */
		for(int kkk = 0; kkk < 10; kkk++){
			matL = matDAD.IC_decomp(diagD, accela_val);
			bool ok = true;
			for(slv_int i = 0; i < size0; i++){
				ok &= (diagD[i] > 0);
			}
			if(ok){
				break;
			}
			accela_val += 0.05;
		}
	} else{
		matL = matDAD.IC_decomp(diagD, accera);
	}
	/* �m�� */
	SparseMat matL_tr = matL.trans();

	/* ���� */
	bool bl = solveICMRTR(size0, conv_cri, max_ite, normB2, diagD, matDAD.matrix, matL.matrix, matL_tr.matrix, vecB2, results, init);

	/* ���ɖ߂� */
	double* result_true = matD*results;
	for(int i = 0; i < size0; i++){
		results[i] = result_true[i];
	}
	delete[] diagD;
	delete[] vecB2;		
	delete[] result_true;
	return bl;

}
/*--------------------*/

/*//=======================================================
// �� MRTR�ŉ���(�E�Ӄm���������v�Z�p�^�[��)
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init){
	const slv_int size = size0;
	dcomplex norm = 0;
	for(int i = 0 ; i < size ; i++){
		norm += vecB[i]*vecB[i];
	}
	double normB2 = abs(norm);
	normB2 = sqrt(normB2);
	return solveICMRTR(size, conv_cri, max_ite, accera, normB2, matA, vecB, results, init);
}

/*//=======================================================
// �� MRTR�ŉ���(���͉E�ӂ�Vector)
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const std::vector<dcomplex>& vecB, std::vector<dcomplex>& results, bool init){
	dcomplex* vecBa = new dcomplex[size0];
	dcomplex* results_a = new dcomplex[size0];
	dcomplex norm = 0;
	for(int i = 0 ; i < size0 ; i++){
		const dcomplex temp = vecB[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = results[i];
	}
	double normB2 = abs(norm);
	normB2 = sqrt(normB2);
	bool bl = solveICMRTR(size0, conv_cri, max_ite, accera, normB2, matA, vecBa, results_a, init);
	for(int i = 0 ; i < size0 ; i++){
		results[i] = results_a[i];
	}
	delete[] results_a;
	delete[] vecBa;
	return bl;
}

/*//=======================================================
// �� MRTR�ŉ���(���͉E�ӂ�Eigen)
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const Eigen::VectorXcd& vecB, dcomplex *results, bool init){
	dcomplex* vecBa = new dcomplex[size0];
	dcomplex norm = 0;
	for(int i = 0 ; i < size0 ; i++){
		const dcomplex temp = vecB(i);
		vecBa[i] = temp;
		norm += temp*temp;
	}
	double normB2 = abs(norm);
	normB2 = sqrt(normB2);
	bool bl = solveICMRTR(size0, conv_cri, max_ite, accera, normB2, matA, vecBa, results, init);
	delete[] vecBa;
	return bl;
}


/*//=======================================================
// �� MRTR�ŉ����E�O�����s�{��
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init){
	/* �Ίp�X�P�[�����O���� */
	if(is_diag_scale){
		return solveICMRTR_diag(size0, conv_cri, max_ite, accera, normB, matA, vecB, results, init);
	}
	/* �R���X�L�[�p�X�p�[�X�s��쐬 */
	dcomplex* diagD = new dcomplex[size0];
	double accela_val = accera;
	/* �����W�������Ȃ玩�����胂�[�h�� */
	SparseMatC matL;
	if(accera < -1){
		accela_val = 1.05;
		/* �Ίp�����ɂȂ�܂Ŏ��{ */
		for(int kkk = 0; kkk < 10; kkk++){
			matL = matA.IC_decomp(diagD, accela_val);
			bool ok = true;
			for(slv_int i = 0; i < size0; i++){
				ok &= (diagD[i].real() > 0);
			}
			if(ok){
				break;
			}
			accela_val += 0.05;
		}
	} else{
		matL = matA.IC_decomp(diagD, accera);
	}
	/* �m�� */
	SparseMatC matL_tr = matL.trans();

	bool bl= solveICMRTR(size0, conv_cri, max_ite, normB, diagD, matA.matrix, matL.matrix, matL_tr.matrix, vecB, results, init);
	delete[] diagD;
	return bl;
}


/*//=======================================================
// �� MRTR�ŉ����E�O�����s�{��
//=======================================================*/
bool MatSolvers::solveICMRTR_diag(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init){
	/* �Ίp�X�P�[�����O */
	dcomplex* vecB2 = new dcomplex[size0];
	SparseMatC matD = matA.diagScaling(vecB2, vecB);
	SparseMatC matDAD = matD*matA*matD;

	dcomplex normBb = 0;
	for(int i = 0; i < size0; i++){
		normBb += vecB2[i]*vecB2[i];
	}
	double normBc = abs(normBb);
	double normB2 = sqrt(normBc);


	/* �R���X�L�[�p�X�p�[�X�s��쐬 */
	dcomplex* diagD = new dcomplex[size0];
	double accela_val = accera;
	/* �����W�������Ȃ玩�����胂�[�h�� */
	SparseMatC matL;
	if(accera < -1){
		accela_val = 1.05;
		/* �Ίp�����ɂȂ�܂Ŏ��{ */
		for(int kkk = 0; kkk < 10; kkk++){
			matL = matDAD.IC_decomp(diagD, accela_val);
			bool ok = true;
			for(slv_int i = 0; i < size0; i++){
				ok &= (diagD[i].real() > 0);
			}
			if(ok){
				break;
			}
			accela_val += 0.05;
		}
	} else{
		matL = matA.IC_decomp(diagD, accera);
	}
	/* �m�� */
	SparseMatC matL_tr = matL.trans();

	bool bl= solveICMRTR(size0, conv_cri, max_ite, normB2, diagD, matDAD.matrix, matL.matrix, matL_tr.matrix, vecB2, results, init);

	/* ���ɖ߂� */
	dcomplex* result_true = matD*results;
	for(int i = 0; i < size0; i++){
		results[i] = result_true[i];
	}
	delete[] diagD;
	delete[] vecB2;		
	delete[] result_true;
	return bl;

}
/*--------------------*/

/*--------------------*/

/*--------------------*/

/*========================================*/
/*========================================*/
/*//=======================================================
// �� MRTR�ŉ����i�{�́j
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
	const double* diagD, const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *vecB, double *results, bool init){

	/* �v�f�m�� */
	Eigen::VectorXd EvecP = Eigen::VectorXd::Zero(size);
	Eigen::VectorXd EvecR(size);
	Eigen::VectorXd EvecU(size);
	Eigen::VectorXd EvecX(size);	
	Eigen::VectorXd EvecY(size);
	Eigen::VectorXd EvecZ(size);


	/* �����ݒ� */
	if(init){
		for(int i = 0 ; i < size ; i++){
			results[i] = 0;
			EvecX(i) = 0.0;
		}
	} else{
		for(int i = 0 ; i < size ; i++){
			EvecX(i) = results[i];
		}
	}

	slv_int* start_posA = new slv_int[size];
	slv_int* end_posA = new slv_int[size];
	matA.getCols(start_posA, end_posA);
	auto col_ptrA = matA.getColPtr();
	auto val_ptrA = matA.getValuePtr();
#ifdef OMP_USING_ICCG		
#pragma omp parallel for
#endif
	for(slv_int ii = 0 ; ii < size ; ii++){
		const slv_int c_size = end_posA[ii];
		double ap_temp=0;
		for(slv_int j = start_posA[ii] ; j < c_size ; j++){
			ap_temp += val_ptrA[j] * results[col_ptrA[j]];
		}
		/* �����c�����v�Z*/
		EvecR(ii) = vecB[ii] - ap_temp;
	}

	/* �ŏ����瓚���������牽�����Ȃ� */
	const double first_normR = EvecR.norm() / normB;
	if(first_normR < conv_cri*0.01){
		delete[] start_posA;
		delete[] end_posA;
		return true;
	}

	/* �O����u=M^-1 * r */
	IC_frbc_process(size, matL, matL_tr, diagD, EvecR, EvecU);	
	/* y0 = -r0 */
	EvecY = -1.0*EvecR;
	/* �O����z=M^-1 * y */
	IC_frbc_process(size, matL, matL_tr, diagD, EvecR, EvecZ);

	/* �ŗǌ��ʂ̕ۑ��p�i�t���O��on�Ȃ�j */
	double* best_results=nullptr;
	double best_resi_value = 1.0e+6;
	if(is_save_best){
		best_results = new double[size];
	}

	bool is_conv = false;
	int bad_counter = 0;

	double zeta=1.0;
	double zeta_old=1.0;
	double eta;
	double nu=1.0;
	/* �����J�n */
	int It = 0;
	for(It = 0; It < max_ite; It++){
		/* v = Au */
		Eigen::VectorXd vecV = matA.matrix * EvecU;
		/* w = M^-1 * v */
		Eigen::VectorXd vecW(size);
		IC_frbc_process(size, matL, matL_tr, diagD, vecV, vecW);

		/* (w, r(k)) */
		double w_r = vecW.dot(EvecR);
		/* (v, w) */
		double v_w = vecV.dot(vecW);

		if(It == 0){
			/* ��(0) */
			zeta = w_r / v_w;
			zeta_old = zeta;
			/* ��(0) */
			eta = 0.0;
		} else{
			/* (w, y(k)) */
			double w_y = vecW.dot(EvecY);
			/* ��(k), ��(k) �̎��̕��� */
			const double temp0 = 1.0 / (nu * v_w - w_y * w_y);
			/* ��(k) */
			zeta = nu * w_r * temp0;
			/* ��(k) */
			eta = -1.0*w_y * w_r * temp0;
		}

		/* ��(k + 1) */
		nu = zeta * w_r;


		/* p(k) */
		double temp = eta * zeta_old / zeta;
		EvecP = EvecU + temp * EvecP;
		zeta_old = zeta;

		/* x(k + 1) */
		EvecX += zeta*EvecP;
		/* y(k + 1) */
		EvecY = eta*EvecY + zeta*vecV;
		/* r(k + 1)*/
		EvecR -= EvecY;

		/* (r(k + 1), r(k + 1)) */
		double norm_r = EvecR.norm();
		/* �������� */
		const double normR = norm_r / normB;
		/* �t���O��on�Ȃ�A�c���ۑ� */
		if(is_save_residual_log){
			residual_log.push_back(normR);
		}
		/* �������� */
		if(normR < conv_cri){
			//std::cout << "Solved!!! -- " << normR  << " " << It << ", " << temp << ", " << temp2  <<std::endl;
			is_conv = true;
			break;
			//}else if( It % 100 == 0 ){
			//	std::cout << "NORMR is -- " << normR  << " " << It << std::endl;
		}
		if(normR < best_resi_value){
			best_resi_value = normR;
			/* �ŗǒl�̍X�V(�t���O��on�Ȃ�) */		
			if(is_save_best){
				for(slv_int i = 0; i < size; i++){
					best_results[i] = results[i];
				}
			}
		}
		/* ���U����P */
		if(diverge_judge_type == 1){
			/* �ŗǒl�~val�ȉ��Ȃ�A���U�J�E���g������ */
			if(normR < best_resi_value * bad_div_val){
				bad_counter = 0;
			}
			/* �ŗǒl�~val�ȏ�Ȃ�A���U�J�E���g�{ */
			if(normR >= best_resi_value * bad_div_val){
				bad_counter++;
			}
			/* ���U�J�E���g��臒l�I�[�o�[�����U�����ŏI��� */
			if(bad_counter >= bad_div_count_thres){
				is_conv = false;
				break;
			}
		}

		/* z(k + 1) */
		EvecZ = eta*EvecZ + zeta*vecW;
		/* u(k + 1) */
		EvecU -= EvecZ;
	}
	delete[] start_posA;
	delete[] end_posA;

	/* ���ʑ�� */
	for(slv_int i = 0; i < size; i++){
		results[i] = EvecX(i);
	}	
	if(!is_conv){
		std::cout << "not Convergence!!! " << std::endl;
		/* �ŗǒl����(�t���O��on�Ȃ�) */
		if(is_save_best){
			for(slv_int i = 0; i < size; i++){
				results[i] = best_results[i];
			}
			delete[] best_results;
		}
	}
	return is_conv ;
}


/*--------------------*/

/*========================================*/
/*========================================*/
/*//=======================================================
// �� CO - MRTR�ŉ����i�{�́j
//=======================================================*/
bool MatSolvers::solveICMRTR(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
	const dcomplex* diagD, const SparseMatBaseC& matA, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *vecB, dcomplex *results, bool init){

	/* �v�f�m�� */
	//double alpha;
	//double beta;
	Eigen::VectorXcd EvecP = Eigen::VectorXcd::Zero(size);
	Eigen::VectorXcd EvecR(size);
	Eigen::VectorXcd EvecU(size);
	Eigen::VectorXcd EvecX(size);	
	Eigen::VectorXcd EvecY(size);
	Eigen::VectorXcd EvecZ(size);


	/* �����ݒ� */
	if(init){
		for(int i = 0 ; i < size ; i++){
			results[i] = 0;
			EvecX(i) = 0.0;
		}
	} else{
		for(int i = 0 ; i < size ; i++){
			EvecX(i) = results[i];
		}
	}

	slv_int* start_posA = new slv_int[size];
	slv_int* end_posA = new slv_int[size];
	matA.getCols(start_posA, end_posA);
	auto col_ptrA = matA.getColPtr();
	auto val_ptrA = matA.getValuePtr();
#ifdef OMP_USING_ICCG		
#pragma omp parallel for
#endif
	for(slv_int ii = 0 ; ii < size ; ii++){
		const slv_int c_size = end_posA[ii];
		dcomplex ap_temp=0;
		for(slv_int j = start_posA[ii] ; j < c_size ; j++){
			ap_temp += val_ptrA[j] * results[col_ptrA[j]];
		}
		/* �����c�����v�Z*/
		EvecR(ii) = vecB[ii] - ap_temp;
	}

	/* �ŏ����瓚���������牽�����Ȃ� */
	const double first_normR = EvecR.norm() / normB;
	if(first_normR < conv_cri*0.01){
		delete[] start_posA;
		delete[] end_posA;
		return true;
	}

	/* �O����u=M^-1 * r */
	IC_frbc_process(size, matL, matL_tr, diagD, EvecR, EvecU);	
	/* y0 = -r0 */
	EvecY = -1.0*EvecR;
	/* �O����z=M^-1 * y */
	IC_frbc_process(size, matL, matL_tr, diagD, EvecR, EvecZ);

	/* �ŗǌ��ʂ̕ۑ��p�i�t���O��on�Ȃ�j */
	dcomplex* best_results=nullptr;
	double best_resi_value = 1.0e+6;
	if(is_save_best){
		best_results = new dcomplex[size];
	}

	bool is_conv = false;
	int bad_counter = 0;

	dcomplex zeta;
	dcomplex zeta_old;
	dcomplex eta;
	dcomplex nu=1.0;
	/* �����J�n */
	int It = 0;
	for(It = 0; It < max_ite; It++){
		/* v = Au */
		Eigen::VectorXcd vecV = matA.matrix * EvecU;
		/* w = M^-1 * v */
		Eigen::VectorXcd vecW(size);
		IC_frbc_process(size, matL, matL_tr, diagD, vecV, vecW);

		/* (w, r(k)) */
		dcomplex w_r = vecW.dot(EvecR);
		/* (v, w) */
		dcomplex v_w = vecV.dot(vecW);

		if(It == 0){
			/* ��(0) */
			zeta = w_r / v_w;
			zeta_old = zeta;
			/* ��(0) */
			eta = 0.0;
		} else{
			/* (w, y(k)) */
			dcomplex w_y = vecW.dot(EvecY);
			/* ��(k), ��(k) �̎��̕��� */
			const dcomplex temp0 = 1.0 / (nu * v_w - w_y * w_y);
			/* ��(k) */
			zeta = nu * w_r * temp0;
			/* ��(k) */
			eta = -1.0*w_y * w_r * temp0;
		}

		/* ��(k + 1) */
		nu = zeta * w_r;


		/* p(k) */
		dcomplex temp = eta * zeta_old / zeta;
		EvecP = EvecU + temp * EvecP;
		zeta_old = zeta;

		/* x(k + 1) */
		EvecX += zeta*EvecP;
		/* y(k + 1) */
		EvecY = eta*EvecY + zeta*vecV;
		/* r(k + 1)*/
		EvecR -= EvecY;

		/* (r(k + 1), r(k + 1)) */
		double norm_r = EvecR.norm();
		/* �������� */
		const double normR = norm_r / normB;
		/* �t���O��on�Ȃ�A�c���ۑ� */
		if(is_save_residual_log){
			residual_log.push_back(normR);
		}
		/* �������� */
		if(normR < conv_cri){
			//std::cout << "Solved!!! -- " << normR  << " " << It << ", " << temp << ", " << temp2  <<std::endl;
			is_conv = true;
			break;
			//}else if( It % 100 == 0 ){
			//	std::cout << "NORMR is -- " << normR  << " " << It << std::endl;
		}
		if(normR < best_resi_value){
			best_resi_value = normR;
			/* �ŗǒl�̍X�V(�t���O��on�Ȃ�) */		
			if(is_save_best){
				for(slv_int i = 0; i < size; i++){
					best_results[i] = results[i];
				}
			}
		}
		/* ���U����P */
		if(diverge_judge_type == 1){
			/* �ŗǒl�~val�ȉ��Ȃ�A���U�J�E���g������ */
			if(normR < best_resi_value * bad_div_val){
				bad_counter = 0;
			}
			/* �ŗǒl�~val�ȏ�Ȃ�A���U�J�E���g�{ */
			if(normR >= best_resi_value * bad_div_val){
				bad_counter++;
			}
			/* ���U�J�E���g��臒l�I�[�o�[�����U�����ŏI��� */
			if(bad_counter >= bad_div_count_thres){
				is_conv = false;
				break;
			}
		}

		/* z(k + 1) */
		EvecZ = eta*EvecZ + zeta*vecW;
		/* u(k + 1) */
		EvecU -= EvecZ;
	}
	delete[] start_posA;
	delete[] end_posA;
	
	/* ���ʑ�� */
	for(slv_int i = 0; i < size; i++){
		results[i] = EvecX(i);
	}	
	if(!is_conv){
		std::cout << "not Convergence!!! " << std::endl;
		/* �ŗǒl����(�t���O��on�Ȃ�) */
		if(is_save_best){
			for(slv_int i = 0; i < size; i++){
				results[i] = best_results[i];
			}
			delete[] best_results;
		}
	}
	return is_conv ;
}




/* end of namespace */
};

