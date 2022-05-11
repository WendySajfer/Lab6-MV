#include <iostream>
#include "windows.h"
#include <vector>
#include "string"
#include <cmath>

using namespace std;

class Matrix {
private:
	vector<vector<double>> matrix;
	vector<int> width_form;
	int m = 0, n = 0;

	void Format() {
		width_form.clear();
		int width, buf_width, width_null;
		string str_width;
		width = 1;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				str_width = to_string(matrix[j][i]);
				for (int i = str_width.size() - 1; i >= 0; i--) {
					if (str_width[i] == '0' || str_width[i] == ',') {
						str_width.erase(i);
					}
					else break;
				}
				buf_width = str_width.size();
				if (width < buf_width) width = buf_width;
			}
			width_form.push_back(width);
		}
	}
	void Size() {
		m = matrix.size();
		n = matrix[0].size();
		for (int i = 1; i < m; i++) {
			if (matrix[i].size() != n) {
				n = 0;
				break;
			}
		}
	}
public:
	void Input_Matrix(int m, int n)
	{
		double number;
		vector<double> str_matrix;
		cout << "Enter the [" << m << "," << n << "] matrix:" << endl;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cin >> number;
				str_matrix.push_back(number);
			}
			matrix.push_back(str_matrix);
			str_matrix.clear();
		}
	}
	void Output_Matrix() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	void Output_SLAU() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
				if (j == n - 2) cout << "| ";
			}
			cout << endl;
		}
		cout << endl;
	}
	vector<vector<double>> get_Matrix() {
		return matrix;
	}
	void set_SLAU(vector<vector<double>> slau) {
		matrix = slau;
	}
};

class Determinant {
private:
	vector<vector<double>> matrix_det;
	double det = 0;
	Matrix M;

	double Accur(double n) {
		double buf = abs(n);
		if (buf != 0) {
			if (fmod(buf, 0.001) < 1) {
				n = round(n * 1000) / 1000;
			}
		}
		return n;
	}
	double Det_2x(vector<vector<double>> matrix_2x) {
		double right = matrix_2x[0][0] * matrix_2x[1][1];
		double left = matrix_2x[0][1] * matrix_2x[1][0];
		return Accur(right - left);
	}
	double Det_3x(vector<vector<double>> matrix_3x) {
		double right = 0;
		double buf = 1, buf1 = 1;
		for (int i = 0; i < 3; i++) {
			buf *= matrix_3x[i][i];
		}
		right += buf;
		buf = 1;
		for (int i = 0; i < 3; i++) {
			if (i + 1 < 3) {
				buf *= matrix_3x[i][i + 1];
				buf1 *= matrix_3x[i + 1][i];
			}
			else {
				buf *= matrix_3x[i][0];
				buf1 *= matrix_3x[0][i];
			}
		}
		right += buf;
		right += buf1;
		double left = 0;
		buf = matrix_3x[0][2] * matrix_3x[2][0] * matrix_3x[1][1];
		left += buf;
		buf = matrix_3x[0][1] * matrix_3x[1][0] * matrix_3x[2][2];
		left += buf;
		buf = matrix_3x[1][2] * matrix_3x[2][1] * matrix_3x[0][0];
		left += buf;
		return Accur(right - left);
	}
	double Det_mx(vector<vector<double>> matrix_mx) {
		if (matrix_mx.size() == 2) return Det_2x(matrix_mx);
		else if (matrix_mx.size() == 3) return Det_3x(matrix_mx);
		else {
			double det_mx = 0;
			for (int i = 0; i < matrix_mx.size(); i++) {
				vector<vector<double>> buf_m = matrix_mx;
				buf_m.erase(buf_m.begin());
				int buf_size = buf_m.size();
				for (int buf_index = 0; buf_index < buf_size; buf_index++) {
					buf_m[buf_index].erase(buf_m[buf_index].begin() + i);
				}
				if (i % 2 == 0) {
					det_mx += matrix_mx[0][i] * Det_mx(buf_m);
				}
				else det_mx -= matrix_mx[0][i] * Det_mx(buf_m);
			}
			return Accur(det_mx);
		}
	}
public:
	void Input_Det_Matrix(vector<vector<double>> matrix) {
		matrix_det = matrix;
	}
	void Decision() {
		det = Det_mx(matrix_det);
		cout << "det(A) = " << det << endl;
	}
	double get_det() {
		return det;
	}
};

class SLAU {
private:
	vector<vector<double>> matrix_slau, matrix_in, A_B_e;
	vector<double> X;
	Matrix M;
	int m;
	void Accur_str(int str_index) {
		double buf;
		vector<double> str_slau = matrix_slau[str_index];
		for (int i = 0; i < str_slau.size(); i++) {
			buf = abs(str_slau[i]);
			if (buf != 0) {
				if (fmod(buf, 0.001) < 1) {
					str_slau[i] = round(str_slau[i] * 1000) / 1000;
				}
			}
		}
		matrix_slau[str_index] = str_slau;
	}
	double Accur(double n) {
		double buf = abs(n);
		if (buf != 0) {
			if (fmod(buf, 0.001) < 1) {
				n = round(n * 1000) / 1000;
			}
		}
		return n;
	}
	void Output_SLAU() {
		M.set_SLAU(matrix_in);
		M.Output_SLAU();
	}
	void X_Record() {
		double buf_x = 0;
		buf_x = Accur((matrix_slau[m - 1][3] - matrix_slau[m - 1][0] * A_B_e[1][m - 2]) / (matrix_slau[m - 1][1] + matrix_slau[m - 1][0] * A_B_e[0][m - 2]));
		X.push_back(buf_x);
		for (int i = m - 2; i >= 0; i--) {
			X.push_back(Accur(A_B_e[0][i] * X[m - i - 2]) + A_B_e[1][i]);
		}
	}
	void Sub_buttom_Str(int str1_index, int str2_index, int null_index) {
		double diff;
		vector<double> str1_slau = matrix_in[str1_index];
		vector<double> str2_slau = matrix_in[str2_index];

		diff = str2_slau[null_index] / str1_slau[null_index];
		for (int i = 0; i < str1_slau.size(); i++) {
			str2_slau[i] -= (str1_slau[i] * diff);
		}
		matrix_in[str2_index] = str2_slau;
	}
	void Sub_top_Str(int str1_index, int str2_index, int null_index) {
		double diff;
		vector<double> str1_slau = matrix_in[str1_index];
		vector<double> str2_slau = matrix_in[str2_index];

		diff = str1_slau[null_index] / str2_slau[null_index];
		for (int i = null_index; i < str2_slau.size(); i++) {
			str1_slau[i] -= (str2_slau[i] * diff);
		}
		matrix_in[str1_index] = str1_slau;
	}
	vector<bool> Test1() {
		bool test1 = true, test2 = true;
		for (int i = 0; i < m; i++) {
			if (matrix_in[i][i] == 0) {
				test1 = false;
				break;
			}
		}
		if (test1) {
			for (int i = 2; i < m; i++) {
				if (matrix_in[0][i] != 0) {
					test2 = false;
					break;
				}
				if (matrix_in[m - 1][i - 2] != 0) {
					test2 = false;
					break;
				}
			}
			for (int i = 1; i < m - 1; i++) {
				for (int j = 0; i < m; i++) {
					if (j != i - 1 && j != i && j != i + 1) {
						if (matrix_in[i][j] != 0) {
							test2 = false;
							break;
						}
					}
				}
			}
		}
		return {test1, test2};
	}
	bool Test2_withCreateSLAU() {
		bool test = true;
		vector<double> str_buf = {0, 0, 0};
		for (int i = 0; i < m; i++) {
			str_buf.push_back(matrix_in[i][m]);
			matrix_slau.push_back(str_buf);
			str_buf.pop_back();
		}
		double a, b, c;
		a = 0;
		for (int i = 0; i < m; i++) {
			if(i != 0) a = matrix_in[i][i - 1];
			b = matrix_in[i][i];
			if (i != m - 1)c = matrix_in[i][i + 1];
			else c = 0;
			if (abs(b) >= (abs(a) + abs(c))) {
				matrix_slau[i][0] = a;
				matrix_slau[i][1] = b;
				matrix_slau[i][2] = c;
			}
			else {
				test = false;
				break;
			}
		}
		return test;
	}
	void Create_3Diagonal() {
		for (int i = 2; i < m; i++) {
			for (int j = 0; j < i - 1; j++) {
				if (matrix_in[i][j] != 0) {
					Sub_buttom_Str(j + 1, i, j);
				}
			}
		}
		for (int i = 0; i < m - 2; i++) {
			for (int j = i + 2; j < m; j++) {
				if (matrix_in[i][j] != 0) {
					Sub_top_Str(i, j - 1, j);
				}
			}
		}
	}
	double ei(int i) {
		return (Accur(matrix_slau[i][0] * A_B_e[0][i - 1] + matrix_slau[i][1]));
	}
	double Ai(int i) {
		return (Accur(matrix_slau[i][0] * A_B_e[0][i - 1] + matrix_slau[i][1]));
	}
	void Create_A_B_e() {
		A_B_e.push_back({ -(Accur(matrix_slau[0][2] / matrix_slau[0][1])) });
		A_B_e.push_back({ Accur(matrix_slau[0][3] / matrix_slau[0][1]) });
		A_B_e.push_back({ 1 });
		for (int i = 1; i < m - 1; i++) {
			A_B_e[2].push_back(ei(i));
			A_B_e[0].push_back(-(Accur(matrix_slau[i][2] / A_B_e[2][i])));
			A_B_e[1].push_back(Accur((matrix_slau[i][3] - matrix_slau[i][0] * A_B_e[1][i - 1]) / A_B_e[2][i]));
		}
	}
public:
	bool Input_SLAU(vector<vector<double>> matrix_A, vector<vector<double>> matrix_B) {
		if (matrix_A.size() != matrix_B.size() || matrix_B[0].size() != 1) {
			return false;
		}
		matrix_in = matrix_A;
		for (int i = 0; i < matrix_A.size(); i++) {
			matrix_in[i].push_back(matrix_B[i][0]);
		}
		m = matrix_A.size();
		vector<bool> t1 = Test1();
		Output_SLAU();
		if (t1[0] == false) return t1[0];
		if (t1[1] == false) {
			Create_3Diagonal();
			Output_SLAU();
			return Test2_withCreateSLAU();
		}
		return Test2_withCreateSLAU();
	}
	void Decision() {
		Create_A_B_e();
		X_Record();
		reverse(X.begin(), X.end());
	}
	void Rezult() {
		cout << "Answer: ";
		for (int i = 0; i < X.size(); i++) {
			cout << "x" << i + 1 << " = " << X[i] << ", ";
		}
		cout << endl;
	}

};

class Task {
private:
	Matrix M;
	Matrix M1;
	SLAU Slau;
	Determinant Det;
	int m, n;
	void Create_SLAU() {
		M.Input_Matrix(m, m);
		Det.Input_Det_Matrix(M.get_Matrix());
		n = 1;
		M1.Input_Matrix(m, n);
	}
public:
	void SLAU_Task() {
		cout << "Input the size of the matrix." << endl << "m = ";
		cin >> m;
		if (m < 3 || m > 49) {
			cout << "Incorrect size of the matrix." << endl;
			return;
		}
		Create_SLAU();
		Det.Decision();
		if (Det.get_det() == 0) {
			cout << "Error! The matrix det = 0." << endl;
			return;
		}
		M.Output_Matrix();
		M1.Output_Matrix();
		bool flag = Slau.Input_SLAU(M.get_Matrix(), M1.get_Matrix());
		if (!flag) {
			cout << "Incorrect matrix." << endl;
			return;
		}
		Slau.Decision();
		Slau.Rezult();
	}
};

int main() {
	setlocale(LC_ALL, "Rus");
	SetConsoleCP(1251);
	int ans, exit = 1;
	while (exit == 1) {
		Task T;
		cout << "1.SLAU" << endl << "2.Exit" << endl << "Choose a way:" << endl;
		cin >> ans;
		switch (ans)
		{
		case 1:
			T.SLAU_Task();
			break;
		case 2:
			exit = 0;
			break;
		default:
			cout << "This task does not exist" << endl;
			break;
		}
	}
	system("pause");
	return 0;
}

/*
3.78 1.08 0
1.08 -2.28 0.37
0 0.37 2.86
0.35
1.27
0.47
*/