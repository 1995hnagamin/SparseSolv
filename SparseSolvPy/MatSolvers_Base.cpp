
#include "MatSolvers.hpp"
#include "SparseMatOperators.hpp"

/* ��p���O��� */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
�\���o��{����
//=======================================================
//=======================================================
*/

/*//=======================================================
// �� �ݒ�R���X�g���N�^
//=======================================================*/
MatSolvers::MatSolvers(){
	is_diag_scale = false;
	is_save_best = false;
	is_save_residual_log = false;
	residual_log.clear();
	/* ���U����F
	0=�ő唽���܂ł��
	1�F�ŗǁ~bad_div_val���傫���l��bad_div_count�����������甭�U�Ƃ��ďI��� */
	diverge_judge_type = 0;
	bad_div_val = 1000.0;
	bad_div_count_thres = 1000;
	/* ��������̐��K���^�C�v(0:�E�ӁA1:�����c���A2:�O���w��) */
	conv_normalize_type = 0;
}

/*//=======================================================
// �� ���O�擾
//=======================================================*/
void MatSolvers::getResidualLog(std::vector<double>& log){
	log.clear();
	const size_t the_size = residual_log.size();
	log.resize(the_size);
	if(is_save_residual_log){
		for(size_t i = 0; i < the_size; i++){
			log[i] = residual_log[i];
		}
	}
}

/*--------------------------------------------------------------------------*/


/*//=======================================================
// �� IC �O�i��ޑ��
//=======================================================*/
void MatSolvers::IC_frbc_process(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const double *vecR, double *vec){
	const slv_int size = size0;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();

	/* ���������v�Z */
	Eigen::VectorXd vec_temp = Eigen::VectorXd(size);
	for (slv_int i = 0; i < size; i++) {
		double s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec_temp(j);
		}
		vec_temp(i) = s / val_ptrL1[c_size];
	}
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
		double s = vec_temp(i);
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] = s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#ifdef AAAAAAAAAAAAA
	const slv_int size = size0;
	//double s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		double s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#endif
}

/*//=======================================================
// �� IC �O�i��ޑ��
//=======================================================*/
void MatSolvers::IC_frbc_process(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();

	Eigen::VectorXd vec_temp = Eigen::VectorXd(size);
	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		double s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec_temp(j);
		}
		vec_temp(i) = s / val_ptrL1[c_size];
	}
	for (slv_int i = size; i-- > 0;) {	
		double s = vec_temp(i);
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) = s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;

#ifdef AAAAAAAAAAAAAAAAAA
	const slv_int size = size0;
	//double s;
	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		double s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#endif
}


/*//=======================================================
// �� IC �O�i��ޑ��
//=======================================================*/
void MatSolvers::IC_frbc_process(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *diagD, const dcomplex *vecR, dcomplex *vec){
	const slv_int size = size0;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	Eigen::VectorXcd vec_temp = Eigen::VectorXcd(size);
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec_temp(j);
		}
		vec_temp(i) = s / val_ptrL1[c_size];
	}
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
		dcomplex s = vec_temp[i];
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] = s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#ifdef AAAAAAAAAAAAAAAAAA
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#endif
}

/*//=======================================================
// �� IC �O�i��ޑ��
//=======================================================*/
void MatSolvers::IC_frbc_process(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *diagD, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	Eigen::VectorXcd vec_temp = Eigen::VectorXcd(size);
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec_temp(j);
		}
		vec_temp(i) = s / val_ptrL1[c_size];
	}
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
		dcomplex s = vec_temp(i);
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) = s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;

#ifdef AAAAAAAAAAAAAAA
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#endif
}


/*--------------------------------------------------------------------------*/



/*//=======================================================
// �� �O�i���
//=======================================================*/
void MatSolvers::fr_process(const slv_int size0, const SparseMatBaseD& matL, const double *vecR, double *vec){
	const slv_int size = size0;
	//double s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		double s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	delete[] start_posL1;
	delete[] end_posL1;
}

/*//=======================================================
// �� �O�i���
//=======================================================*/
void MatSolvers::fr_process(const slv_int size0, const SparseMatBaseD& matL, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;
	//double s;
	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();
	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		double s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	delete[] start_posL1;
	delete[] end_posL1;
}


/*//=======================================================
// �� �O�i���
//=======================================================*/
void MatSolvers::fr_process(const slv_int size0, const SparseMatBaseC& matL, const dcomplex *vecR, dcomplex *vec){
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	delete[] start_posL1;
	delete[] end_posL1;
}

/*//=======================================================
// �� �O�i���
//=======================================================*/
void MatSolvers::fr_process(const slv_int size0, const SparseMatBaseC& matL, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	/* ���������v�Z */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	delete[] start_posL1;
	delete[] end_posL1;
}


/*----------------------------------------------*/



/*//=======================================================
// �� ��ޑ��
//=======================================================*/
void MatSolvers::bc_process(const slv_int size0, const SparseMatBaseD& matL_tr, const double *vecR, double *vec){
	const slv_int size = size0;
	//double s;

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();

	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = vecR[i];
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec[the_row];
		}
		vec[i] = s  / val_ptrL2[start_posL2[i]];
	}
	delete[] start_posL2;
	delete[] end_posL2;
}

/*//=======================================================
// �� ��ޑ��
//=======================================================*/
void MatSolvers::bc_process(const slv_int size0, const SparseMatBaseD& matL_tr, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;
	//double s;

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = EvecR(i);
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec(the_row);
		}
		vec(i) = s  / val_ptrL2[start_posL2[i]];
	}
	delete[] start_posL2;
	delete[] end_posL2;
}


/*//=======================================================
// �� ��ޑ��
//=======================================================*/
void MatSolvers::bc_process(const slv_int size0, const SparseMatBaseC& matL_tr, const dcomplex *vecR, dcomplex *vec){
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = vecR[i];
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec[the_row];
		}
		vec[i] = s  / val_ptrL2[start_posL2[i]];
	}
	delete[] start_posL2;
	delete[] end_posL2;
}

/*//=======================================================
// �� ��ޑ��
//=======================================================*/
void MatSolvers::bc_process(const slv_int size0, const SparseMatBaseC& matL_tr, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* ���������v�Z */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = EvecR(i);
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec(the_row);
		}
		vec(i) = s  / val_ptrL2[start_posL2[i]];
	}
	delete[] start_posL2;
	delete[] end_posL2;
}

/*//=======================================================
// �� �����W���̎�������
//=======================================================*/
void MatSolvers::auto_accel_determine(const slv_int size0, double accel_ini, const SparseMat& matA, double* diagD, SparseMat& matL){
	double accela_val = fabs(accel_ini);
	if(accela_val < 0.9 || accela_val > 1.8){
		accela_val = 1;
	}
	/* �Ίp�����ɂȂ�܂Ŏ��{ */
	for(int kkk = 0; kkk < 80; kkk++){
		SparseMat matL_temp = matA.IC_decomp(diagD, accela_val);
		bool ok = true;
		for(slv_int i = 0; i < size0; i++){
			ok &= (diagD[i] > 0);
		}
		if(ok){
			matL = std::move(matL_temp);
			return;
		}
		accela_val += 0.01;
	}
}

/*//=======================================================
// �� �����W���̎�������(���f)
//=======================================================*/
void MatSolvers::auto_accel_determine(const slv_int size0, double accel_ini, const SparseMatC& matA, dcomplex* diagD, SparseMatC& matL){
	double accela_val = fabs(accel_ini);
	if(accela_val < 0.9 || accela_val > 1.8){
		accela_val = 1;
	}
	/* �Ίp�����ɂȂ�܂Ŏ��{ */
	for(int kkk = 0; kkk < 80; kkk++){
		SparseMatC matL_temp = matA.IC_decomp(diagD, accela_val);
		bool ok = true;
		for(slv_int i = 0; i < size0; i++){
			ok &= (diagD[i].real() > 0);
		}
		if(ok){
			matL = std::move(matL_temp);
			return;
		}
		accela_val += 0.01;
	}
}

/* end of namespace */
};

