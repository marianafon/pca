#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "jama_eig.h"

using namespace std;
using namespace TNT;
using namespace JAMA;

namespace PCA {
	bool debug = false;
	
	template <class T>
	
	void convert_from_string(T& value, const string& s)
	{
		stringstream ss(s);
		ss >> value;
	}

	void load_data_from_file(Array2D<double>& d, char*& file_path) {
		ifstream in(file_path);
        string line;
        int r = 0;
		
        if (in.is_open()) {
            while (in.good()) {
				int col = 0;
                getline(in, line);
                if (line.empty()) continue;
                
				size_t start_pos = 0;
				char space = ','; // modo csv
				while (true) {
					size_t pos = line.find(space, start_pos);
					string data(line.substr(start_pos, pos - start_pos));
					if (!data.empty()) {
						double v = 0;
						convert_from_string(v, data);
						
						if (debug)
							cout << "value: " << v << endl;
						d[r][col] = v;
					}
              
					if ((int)pos != -1) {
						start_pos = pos + 1;
					}
					else {
						break;
					}
					col += 1;
				}
				r += 1;
            }
            in.close();
        }
	}
	
	void adjust_data(Array2D<double>& d, Array1D<double>& means) {
		double mean;
		int i, j;
		for (i=0; i<d.dim2(); ++i) { 
			mean = 0;
			for (j=0; j<d.dim1(); ++j) {
			   mean += d[j][i];
			}
			// calculate the mean itself
			mean /= d.dim1();

			// store the mean per dimension
			means[i] = mean;		   

			// subtract the mean
			for (j=0; j<d.dim1(); ++j) {
			   d[j][i] -= mean;
			}
		}
    }

	double compute_covariance(const Array2D<double>& d, int i, int j) {
		double cov = 0;
		int k;
		// note que 'd' já está reduzido pela média
		for (k=0; k<d.dim1(); ++k) {
			cov += d[k][i] * d[k][j];
		}
		// (sum / n-1)
		return cov / (d.dim1() - 1);
	}

	void compute_covariance_matrix(const Array2D<double> & d, Array2D<double> & covar_matrix) {
		int dim = d.dim2(); // duas dimensões -> gerar matriz de covariância 2x2
        int i, j;
		assert(dim == covar_matrix.dim1());
		assert(dim == covar_matrix.dim2());
        for (i=0; i<dim; ++i) {
			for (j=i; j<dim; ++j) {
				// calcular covariância entre i e j
				covar_matrix[i][j] = compute_covariance(d, i, j);
			}
		}
		// fill the Left triangular matrix
		for (i=1; i<dim; i++) {
			for (j=0; j<i; ++j) {
				covar_matrix[i][j] = covar_matrix[j][i];
			}
		}

	}

	// Calculate the eigenvectors and eigenvalues of the covariance matrix
	void eigen(const Array2D<double> & covar_matrix, Array2D<double>& eigenvector, Array2D<double>& eigenvalue) {
		Eigenvalue<double> eig(covar_matrix);
		eig.getV(eigenvector);
		eig.getD(eigenvalue);
	}


	void transpose(const Array2D<double>& src, Array2D<double>& target) {
        int i,j;
        for (i=0; i<src.dim1(); ++i) {
			for (j=0; j<src.dim2(); ++j) {
				target[j][i] = src[i][j];
			}
		}
	}

	// z = x * y
	void multiply(const Array2D<double>& x, const Array2D<double>& y, Array2D<double>& z) {
		assert(x.dim2() == y.dim1());
        int i, j, k, d;
        double sum = 0;
        for (i=0; i<x.dim1(); ++i) {
			for (j=0; j<y.dim2(); ++j) {
				sum = 0;
				d = y.dim1();
				for (k=0; k<d; k++) {
					sum += x[i][k] * y[k][j];
				}
				z[i][j] = sum;
			}
		}
	}
	
	void choose_components(Array2D<double>& a_valor, Array2D<double>& a_vetor, Array2D<double>& featureVector, int K) {
		// ordenar autovetores a partir dos autovalores.
		// odd-even.
		int n = a_valor.dim1(); // col
		int phase, i, j;
		double temp, temp2;
        for (phase = 0; phase < n; phase++) {
			if (phase % 2 == 0) { // fase par
				for (i = 1; i < n; i += 2) {
					if (a_valor[i-1][i-1] < a_valor[i][i]) {
						temp = a_valor[i][i];
						a_valor[i][i] = a_valor[i-1][i-1];
						a_valor[i-1][i-1] = temp;
						// ordenamos o autovetor equivalente
						for (j = 0; j < n; j++) {
							temp2 = a_vetor[j][i];
							a_vetor[j][i] = a_vetor[j][i-1];
							a_vetor[j][i-1] = temp2;
						} // fim for j
				    } // fim if swap
				} // fim for i
			} // fim if fase par
			else { // fase impar
				for (i = 1; i < n-1; i += 2) {
					if (a_valor[i][i] < a_valor[i+1][i+1]) {
						temp = a_valor[i][i];
						a_valor[i][i] = a_valor[i+1][i+1];
						a_valor[i+1][i+1] = temp;
						// ordenamos o autovetor equivalente
						for (j = 0; j < n; j++) {
							temp2 = a_vetor[j][i];
							a_vetor[j][i] = a_vetor[j][i+1];
							a_vetor[j][i+1] = temp2;
						} // fim for j
					} // fim if swap
				} // fim for i
			} // fim if fase impar
		} // fim for fases
		
		// selecionar os K primeiros autovetores, mantendo o numero de linhas original.
		for (i = 0; i < n; i++) {
			for (j = 0; j < K; j++) {
				featureVector[i][j] = a_vetor[i][j];
			}
		}		
		// resultado final: featureVector.		
	} // fim choose_components function
	
	void normalize(Array2D<double>& data, const Array1D<double>& means) {
        int i, j;
        for (i=0; i<data.dim2(); i++) {
			for (j=0; j<data.dim1(); j++) {
				data[j][i] += means[i];
			}
		}
	}
	
	void print_array(const Array2D<double>& x) {
		for (int i=0; i<x.dim1(); i++) {
			for (int j=0; j<x.dim2(); j++) {
				cout << x[i][j] << "  ";
			}
			cout << endl;
		}
		cout << "------" << endl;
	}
	
	void export_array(const Array2D<double>& x, string filename) {
		// exportar no modo csv
		ofstream outp (filename.c_str());
		for (int i=0; i<x.dim1(); i++) {
			for (int j=0; j<x.dim2(); j++) {
				outp << x[i][j];
				// enquanto não chega no último...
				if (j != x.dim2()-1) { 
					outp << ",";
				}
			}
			outp << endl;
		}
		outp.close();
	}
	
	int get_dimension_from_filename(string filename) {
		int i;
		string dim;
		for (i=filename.size()-5; i>0; i--) {
			if (filename[i] == '/')
				break;
			dim += filename[i];
		}
		reverse(dim.begin(), dim.end());
		int dimension = strtol(dim.c_str(), NULL, 10);
		return dimension;
	}
	
}

int main(int argc, char* argv[]) {
    if (argc != 3) {        
        cout << "Usage: " << argv[0] << " <data_file> <K>" << endl;        
        return -1;   
    }
	using namespace PCA;
	
	string filename = argv[1];
	
	// informação do banco de dados
    int dim = get_dimension_from_filename(filename);
    int K = strtol(argv[2], NULL, 10);
    
    // verificar argumentos corretos	
    assert(dim >= K);
    const int row = 125; // numero de dados por dimensão
    const int col = dim; // numero de dimensões
    //const int col = 2; // numero de dimensões

    Array2D<double> d(row, col); // matriz de dados 
    load_data_from_file(d, argv[1]); // carrega a matriz de dados

    Array1D<double> means(col);
    adjust_data(d, means); // subtrai as médias

    Array2D<double> covar_matrix(col, col);
    compute_covariance_matrix(d, covar_matrix); // computa matriz de covariância
	
    // obter os autovetores e autovalores através da matriz de covariância.
    Array2D<double> eigenvector(col, col);
    Array2D<double> eigenvalue(col, col);
    eigen(covar_matrix, eigenvector, eigenvalue);

	
	// determinar o featureVector reduzindo os eigenvectors principais em K eigenvectors.
	// resultado possui K colunas (dimensões)
	Array2D<double> featureVector(col, K);
	choose_components(eigenvalue, eigenvector, featureVector, K);	

    
    // final_data^T = FeatureVector^T * DataAdjustada^T
	// sendo final_data a matriz calculada com os melhores eigenvectors
    Array2D<double> fd_transp(K, row);
	Array2D<double> final_data(row, K);
	Array2D<double> transpose_featV(K, col);
    Array2D<double> transpose_data(col, row);
	
	transpose(featureVector, transpose_featV);
    transpose(d, transpose_data);
    multiply(transpose_featV, transpose_data, fd_transp);
	transpose(fd_transp, final_data); // forma original, cada coluna é uma dimensão.	
	
	export_array(final_data, "teste.out");
	// caso queira reverter os dados para as dimensões originais...
	/*
	// gerando dados novos
	Array2D<double> dados_novos_transp(col, row);
	Array2D<double> dados_novos(row, col);
	
	// dados_novos^T = (featureVector * final_data^T) + means
	multiply(featureVector, fd_transp, dados_novos_transp);
	transpose(dados_novos_transp, dados_novos);
	// somar com as médias, pra retornar ao valor inicial
	normalize(dados_novos, means);
    
	ostringstream filepath;
	string s;
	filepath << "./outputs/" << dim << "to" << K << ".out";
	s = filepath.str();
	export_array(dados_novos, s);
	*/
    return 0;
}
