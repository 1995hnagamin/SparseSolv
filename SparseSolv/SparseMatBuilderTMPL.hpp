#ifndef DEF_SPR_MAT_BUILD_TMPL
#define DEF_SPR_MAT_BUILD_TMPL

#include "SparseMatTMPL.hpp"

namespace SRLfem {

template<typename DType>
class SparseMatBuilderTMPL {
private:
	slv_int size;																/* 行列の行数 */
	std::unique_ptr<std::map<slv_int, DType>[]> tempMat;						/* 一時保存行列（キーが列位置） */
public:
	SparseMatBuilderTMPL();
	SparseMatBuilderTMPL(slv_int size0);										/* コンストラクタ */
	slv_int getSize() const{ return size; };

	void tempInitialize(slv_int ss);											/* 一時行列を作成 */
	void tempInitialize();														/* 一時行列を作成 */
	void resizeInitialize(slv_int new_size);									/* 一時行列を作成&初期化 */

	SparseMatTMPL<DType> build(bool toSquare=false);							/* 一時配列を確定させる(bool toSquare: trueなら0を入れて正方行列にする) */
	bool isInclude(slv_int gyo, slv_int target_r, slv_int& result_retu) const;	/* i行目にtarget_r列があるかどうか（あったらそのindexを返す） */	
	slv_int getMaxCol() const;													/* スパース内の最大の列位置を返す */
	void add(slv_int gyo, slv_int retu, DType val);								/* 一時配列にpush */

	void readMat(const std::string& filename);
};

/*=======================================================================*/
/*=======================================================================*/

/*//=======================================================
// ● コンストラクタ
//=======================================================*/
template<typename DType>
SparseMatBuilderTMPL<DType>::SparseMatBuilderTMPL():
		size(0), tempMat(nullptr){
}

template<typename DType>
SparseMatBuilderTMPL<DType>::SparseMatBuilderTMPL(slv_int size0):
		size(size0), tempMat(nullptr){
	/* 一時行列初期化 */
	tempInitialize();
}


/*//=======================================================
// ● 一時行列を作成
//=======================================================*/
template<typename DType>
void SparseMatBuilderTMPL<DType>::tempInitialize(slv_int ss){
	size = ss;
	this->tempInitialize();
}

/*//=======================================================
// ● 一時行列の初期化
//=======================================================*/
template<typename DType>
void SparseMatBuilderTMPL<DType>::tempInitialize(){
	assert(size > 0 && "Mat size is not initialized ! ");
	tempMat = std::make_unique<std::map<slv_int, DType>[]>(size);
	for(slv_int i = 0 ; i < size ; i++){
		tempMat[i].clear();
	}
}

/*//=======================================================
// ● 一時行列を作成&初期化
//=======================================================*/
template<typename DType>
void SparseMatBuilderTMPL<DType>::resizeInitialize(slv_int new_size){
	tempMat.reset(nullptr);
	size = new_size;
	tempMat = std::make_unique<std::map<slv_int, DType>[]>(size);
	for(slv_int i = 0 ; i < size ; i++){
		tempMat[i].clear();
	}
}

/*//=======================================================
// ● 一時配列を確定させる
//=======================================================*/
template<typename DType>
SparseMatTMPL<DType> SparseMatBuilderTMPL<DType>::build(bool toSquare){
	/* Eigenに代入 */
	slv_int max_retu = 0;
	std::vector<Eigen::Triplet<DType>> tripletList;
	for(slv_int i = 0; i < size; i++) {
		const slv_int ss = tempMat[i].size();
		if(ss == 0) {
			continue;
		}
		for(auto& itr : tempMat[i]) {
			const slv_int pos = itr.first;
			const DType val = itr.second;
			tripletList.push_back(Eigen::Triplet<DType>(i, pos, val));
			if(pos > max_retu){
				max_retu = pos;
			}
		}
		tempMat[i].clear();
	}
	max_retu++;
	/* Eigenに確定(列方向の最大値も指定する必要あり)。 */
	/* toSquare=true なら、行数の方が大きいなら正方サイズにする */
	/* toSquare=false なら、長方形のまま */
	if(toSquare){		
		if(max_retu < size){
			tripletList.push_back(Eigen::Triplet<DType>(size-1, size-1, 0.0));
			max_retu = size;
		}	
	}
	Eigen::SparseMatrix<DType, Eigen::RowMajor> matrix(size, max_retu);
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
	return matrix;
}

/*//=======================================================
// ● i行目にtarget_r列があるかどうか（あったらそのindexを返す）
//=======================================================*/
template<typename DType>
bool SparseMatBuilderTMPL<DType>::isInclude(slv_int gyo, slv_int target_r, slv_int& result_retu) const{
	slv_int counter=0;
	for(const auto& itr : tempMat[gyo]){
		if( itr.first == target_r){				
			result_retu = counter;
			return true;
		}
		counter++;
	}
	result_retu = -1;
	return false;
}

/*//=======================================================
// ● mapにインサート
//=======================================================*/
template<typename DType>
void SparseMatBuilderTMPL<DType>::add(slv_int gyo, slv_int retu, DType val){
	tempMat[gyo][retu] += val;
}

 
/*//=======================================================
// ● スパース内の最大の列位置を返す
//=======================================================*/
template<typename DType>
slv_int SparseMatBuilderTMPL<DType>::getMaxCol()const{
	slv_int max = 0;
	for(slv_int i = 0; i < size; i++) {
		if(tempMat[i].empty()) {
			continue;
		}
		const slv_int tmp = tempMat[i].crbegin()->first;
		if(max < tmp){
			max = tmp;
		}
	}
	return max;
}


/*//=======================================================
// ● ファイルから行列構成
//=======================================================*/
template<typename DType>
void SparseMatBuilderTMPL<DType>::readMat(const std::string& filename){
	std::fstream fp(filename, std::ios::in);
	slv_int read_size;
	std::string temp_str1, temp_str2, temp_str3;
	fp >> temp_str1 >> temp_str2  >> temp_str3 >> read_size;
	this->size = read_size;
	std::cout << "read size " << read_size << std::endl;
	tempMat = std::make_unique<std::map<slv_int, DType>[]>(size);

	int temp1, temp2;
	slv_int* temp_sizes = new slv_int[size];
	slv_int** temp_cols = new slv_int*[size];
	for(slv_int i = 0 ; i < size ; i++){
		fp >> temp_sizes[i];
		temp_cols[i] = new slv_int[temp_sizes[i]];
		std::cout << temp_sizes[i] << std::endl;
	}
	/**/
	/**/
	fp >> temp_str1 >> temp_str2;
	for(slv_int i = 0 ; i < size ; i++){
		if(temp_sizes[i] == 0){
			fp >> temp1;
		}else{
			for(slv_int j = 0; j < temp_sizes[i]; j++){
				fp >> temp2;
				temp_cols[i][j] = temp2;
				std::cout << i << ", " << j << " >" << temp2 << std::endl;
			}
		}
	}
	/**/
	/**/
	fp >> temp_str1 >> temp_str2  >> temp_str3;
	std::cout << temp_str1 << ", " << temp_str2  << ", " <<  temp_str3 << std::endl;;
	for(slv_int i = 0 ; i < size ; i++){
		for(slv_int j = 0; j < temp_sizes[i]; j++){
			double tempD;
			fp >> tempD;
			std::cout << i << ", " << j << " >" << tempD << std::endl;
			this->add(i, temp_cols[i][j], tempD);
		}
	}
	fp.close();
	//
	delete[] temp_sizes;	
	for(slv_int i = 0 ; i < size ; i++){
		delete[] temp_cols[i];
	}
	delete[] temp_cols;
}

} /* namespace SRLfem */

#endif
