/***********************************************************/
/**	\file
	\brief		Read and Write of Math Related Materials
	\details	matrix files are stored 1-indexing
	\author		Yizhong Zhang
	\date		9/19/2012
*/
/***********************************************************/
#ifndef __YZ_MATH_IO_H__
#define __YZ_MATH_IO_H__

#include <iostream>
#include <fstream>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"

namespace yz{

//	========================================
///@{
/**	@name Read & Write Sparse Matrix
*/
//	========================================

/**
	Read sparse matrix from 1-indexing file

	\param	matrix_file_name	file name of the matrix file
	\param	spa_mat				the matrix to store the reading data
	\return						whether read succeed
*/
template<typename T>
int readSparseMatrixFromFile(const char*			matrix_file_name,
							 SparseMatrixBase<T>&	spa_mat){
	std::ifstream file;
	file.open(matrix_file_name);

	if( !file.is_open() ){
		#ifndef BE_QUIET
			std::cout << "error: readSparseMatrixFromFile, cannot open " << matrix_file_name << std::endl;
		#endif
		return 0;
	}

	spa_mat.Reset();
	int row, col, nnz;
	T val;
	file >> row >> col >> nnz;

	spa_mat.SetDimension(row, col);
	while( nnz-- ){
		file >> row >> col >> val;
		row --;
		col --;
		spa_mat.InsertElement(val, row, col);
	}

	spa_mat.Reorder();

	file.close();
	return 1;
}

/**
	Write sparse matrix to 1-indexing file

	symmetric matrix is stored fully

	\param	matrix_file_name	file name of the matrix file
	\param	spa_mat				the matrix
	\return						whether write succeed
*/
template<typename T>
int writeSparseMatrixToFile(const char*			matrix_file_name,
							SparseMatrixCoo<T>&	spa_mat){
	std::ofstream file;
	file.open(matrix_file_name);

	if( !file.is_open() ){
		#ifndef BE_QUIET
			std::cout << "error: WriteSparseMatrixToFile, cannot open " << matrix_file_name << std::endl;
		#endif
		return 0;
	}

	file << spa_mat.row_num << '\t' << spa_mat.col_num << '\t' << spa_mat.NNZ() << std::endl;

	for(int i=0; i<spa_mat.NNZ(); i++){
		file << spa_mat.row_id[i]+1 << '\t' << spa_mat.col_id[i]+1 << '\t' << spa_mat.value[i] << std::endl;
		if( spa_mat.IsSymmetric() && spa_mat.row_id[i] != spa_mat.col_id[i] )	//	non diagonal elements, switch row and col and output
			file << spa_mat.col_id[i]+1 << '\t' << spa_mat.row_id[i]+1 << '\t' << spa_mat.value[i] << std::endl;
	}

	file.close();
	return 1;
}

/**
	Write sparse matrix to 1-indexing file

	symmetric matrix is stored fully

	\param	matrix_file_name	file name of the matrix file
	\param	spa_mat				the matrix
	\return						whether write succeed
*/
template<typename T>
int writeSparseMatrixToFile(const char*			matrix_file_name,
							SparseMatrixCSR<T>&	spa_mat){
	std::ofstream file;
	file.open(matrix_file_name);

	if( !file.is_open() ){
		#ifndef BE_QUIET
			std::cout << "error: writeSparseMatrixToFile, cannot open " << matrix_file_name << std::endl;
		#endif
		return 0;
	}

	file << spa_mat.row_num << '\t' << spa_mat.col_num << '\t' << spa_mat.NNZ() << std::endl;

	for(int i=0; i<spa_mat.row_num; i++){
		for(int k=spa_mat.row_start[i]; k<spa_mat.row_start[i+1]; k++){
			int j = spa_mat.col_id[k];
			file << i+1 << '\t' << j+1 << '\t' << spa_mat.value[k] << std::endl;
			if( spa_mat.IsSymmetric() && i!=j )
				file << j+1 << '\t' << i+1 << '\t' << spa_mat.value[k] << std::endl;
		}
	}

	file.close();
	return 1;
}

///@}

//	========================================
///@{
/**	@name Read & Write Dense Vector
*/
//	========================================

/**
	Read Dense Vector from file

	\param	vector_file_name	file name of the vector file
	\param	dense_vec			the vector to store the reading data
	\return						whether read succeed
*/
template<typename T>
int readDenseVectorFromFile(const char*		vector_file_name,
							DenseVector<T>&	dense_vec){
	std::ifstream file;
	file.open(vector_file_name);

	if( !file.is_open() ){
		#ifndef BE_QUIET
			std::cout << "error: readDenseVectorFromFile, cannot open " << vector_file_name << std::endl;
		#endif
		return 0;
	}

	dense_vec.Reset();

	TYPE_PROMOTE(T, float) val;
	while( file >> val ){
		dense_vec.value.push_back( val );
	}

	file.close();
	return 1;
}

/**
	Write Dense Vector to File

	\param	vector_file_name	file name of the vector file
	\param	dense_vec			the vector 
	\return						whether write succeed
*/
template<typename T>
int writeDenseVectorToFile(const char*		vector_file_name,
						   DenseVector<T>&	dense_vec){
	std::ofstream file;
	file.open(vector_file_name);

	if( !file.is_open() ){
		#ifndef BE_QUIET
			std::cout << "error: writeDenseVectorToFile, cannot open " << vector_file_name << std::endl;
		#endif
		return 0;
	}

	for(int i=0; i<dense_vec.Dim(); i++){
		file << dense_vec[i] << '\t';
	}

	file.close();
	return 1;
}



}	//	end namespace yz

#endif	//	__YZ_MATH_IO_H__