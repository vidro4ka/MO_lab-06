#include<vector>
#include<iostream>
#include<string>
#include<math.h>
#include<iomanip>
#include<map>
#include<exception>

double null(double number) {
	if (number == -0) {
		number = 0;
		return number;
	}
	return number;
}

int nakop = 1;

void printer(std::vector<std::vector<double>>& matr, const std::vector<std::string>& basis,
	const std::vector<std::string>& free) {
	std::cout << "\n\n\t ";
	for (size_t i = 0; i < free.size(); ++i) {
		std::cout << free[i] << "	 ";
	}
	std::cout << "\n";
	for (size_t i = 0; i < basis.size(); ++i) {
		std::cout << basis[i] << "	";
		for (size_t j = 0; j < matr[i].size(); ++j) {
			if (matr[i][j] == -0) {
				matr[i][j] = 0;
			}
			std::cout << std::setw(6) << std::fixed << std::setprecision(3) << round(matr[i][j] * 1000) / 1000 << std::right << "\t";
		}
		std::cout << "\n";
	}
}

int find_column(std::vector<double>& F) {
	int index_r_column = -1;
	for (size_t j = 1; j < F.size(); ++j) {
		if (F[j] > 0) {
			index_r_column = j;
		}
	}
	return index_r_column;
}

int find_row(std::vector<std::vector<double>>& sympl, int index_r_column) {
	int index_r_row = -1;
	double min = INT64_MAX;
	for (size_t i = 0; i < sympl.size() - 1; ++i) {
		if (sympl[i][0] >= 0 && sympl[i][index_r_column] < 0) {
			continue;
		}
		double znach = sympl[i][0] / sympl[i][index_r_column];
		if (znach >= 0 && znach < min) {
			min = znach;
			index_r_row = i;
		}
	}

	return index_r_row;
}

void transformation(std::vector<std::vector<double>>& sympl, int index_r_row, int index_r_column) {
	double razr = sympl[index_r_row][index_r_column];
	for (size_t i = 0; i < sympl.size(); ++i) {
		if (i != index_r_row) {
			for (size_t j = 0; j < sympl[i].size(); ++j) {
				if (j != index_r_column) {
					sympl[i][j] = sympl[i][j] - (sympl[index_r_row][j] * sympl[i][index_r_column] / razr);
				}
			}
		}
	}
	for (size_t i = 0; i < sympl.size(); ++i) {
		if (i != index_r_row) {
			sympl[i][index_r_column] = -1 * sympl[i][index_r_column] / razr;
		}
	}
	for (size_t j = 0; j < sympl[index_r_row].size(); ++j) {
		if (j != index_r_column) {
			sympl[index_r_row][j] = sympl[index_r_row][j] / razr;
		}
	}
	sympl[index_r_row][index_r_column] = 1 / razr;
}

void searher(std::vector<std::vector<double>>& matr, std::vector<std::string>& free, std::vector<std::string>& basis) {
	int index_r_column;
	for (size_t i = 0; i < matr.size() - 1; ++i) {
		if (matr[i][0] < 0) {
			bool otr = false;
			size_t index_r_raw = i;
			for (size_t j = 1; j < matr[i].size(); ++j) {
				if (matr[i][j] < 0) {
					index_r_column = j;
					//size_t index_r_raw = find_row(matr, index_r_column);
					
					transformation(matr, index_r_raw, index_r_column);
					std::swap(free[index_r_column], basis[index_r_raw]);
					
					otr = true;
					break;
				}
			}
			if (otr) {
				printer(matr, basis, free);
				std::cout << "\n ITERATION\n";
				i = 0;
			}
			else {
				throw std::logic_error("No sollution to this task :(");
			}
		}
	}
}

void method(std::vector<std::vector<double>>& matr, std::vector<std::string>& free,
	std::vector<std::string>& basis, std::string rezult) {

	if (rezult == "min") {
		for (auto& el : matr[matr.size() - 1]) {
			el = -el;
		}
	}

	std::map<std::string, double> m_result;

	for (size_t i = 1; i < free.size(); ++i) {
		m_result[free[i]] = 0;
	}

	printer(matr, basis, free);
	searher(matr, free, basis);

	int index_r_column = find_column(matr[matr.size() - 1]);

	while (index_r_column != -1) {
		int index_r_row = find_row(matr, index_r_column);
		if (index_r_row == -1) {
			std::cout << R"(Infinity number of sollution \(-_-)/)";
			break;
		}
		std::swap(basis[index_r_row], free[index_r_column]);
		transformation(matr, index_r_row, index_r_column);
		printer(matr, basis, free);
		index_r_column = find_column(matr[matr.size() - 1]);
	}
	std::cout << "Answer: " << std::endl;
	for (size_t i = 0; i < basis.size() - 1; ++i) {
		m_result.at(free[i]) = matr[i][0];
	}
	for (const auto& el : m_result) {
		std::cout << el.first << " = " << el.second << std::endl;
	}
	if (rezult == "max") {
		matr[matr.size() - 1][0] *= -1;
	}
	std::cout << "F = " << matr[matr.size() - 1][0] << std::endl;
}

class Simplex_tabels {
	std::vector<std::string> basis;
	std::vector<std::string> free;
	std::vector<std::vector<double>> matr;
	std::string rezult;
	std::map<std::string, double> m_result;
	std::map<std::string, double> the_int_result;
	double promezutochniy_result_func;
	int maximum = INT16_MIN;

public:

	Simplex_tabels(std::vector<std::string>& _basis, std::vector<std::string>& _free,
		std::vector<std::vector<double>>& _matr, std::string _rezult) {
		basis = _basis;
		free = _free;
		matr = _matr;
		rezult = _rezult;
	}

	void Simpl_method() {

		if (rezult == "min") {
			for (auto& el : matr[matr.size() - 1]) {
				el = -el;
			}
		}

		for (size_t i = 1; i < free.size(); ++i) {
			m_result[free[i]] = 0;
		}

		printer(matr, basis, free);
		searher(matr, free, basis);

		int index_r_column = find_column(matr[matr.size() - 1]);

		while (index_r_column != -1) {

			int index_r_row = find_row(matr, index_r_column);

			if (index_r_row == -1) {
				throw std::logic_error("Infinity number of sollution 204 (-_-)/");
			}

			std::swap(basis[index_r_row], free[index_r_column]);

			transformation(matr, index_r_row, index_r_column);

			printer(matr, basis, free);

			index_r_column = find_column(matr[matr.size() - 1]);
		}
		std::cout << "Answer:\n";
		for (size_t i = 0; i < basis.size() - 1; ++i) {
			try {
				m_result.at(basis[i]) = matr[i][0];
			}
			catch (std::exception&) {}
		}
		for (const auto& el : m_result) {
			std::cout << el.first << " = " << el.second << std::endl;
		}
		promezutochniy_result_func = matr[matr.size() - 1][0];
		if (rezult == "min") {
			promezutochniy_result_func *= -1;
		}
		std::cout << "F  = " << promezutochniy_result_func << std::endl;
	}

	Simplex_tabels& operator=(const Simplex_tabels& ob) {
		matr = ob.matr;
		free = ob.free;
		basis = ob.basis;
		rezult = ob.rezult;
		return *this;
	}
};




int main() {
	std::vector<std::vector<double>> matr = {
		{-1, -1, -2, -7},
		{-1, -3, -6, -2},
		{-1, -9, -2, -6},
		{-1, -6, -3, -5},
		{0, -1, -1, -1}
	};
	std::vector<std::string> free = { "sv", "x1", "x2", "x3" };
	std::vector<std::string> basis = { "x4", "x5", "x6", "x7", "F"};
	std::string rezult = "max";

	Simplex_tabels ob(basis, free, matr, rezult);

	try {
		ob.Simpl_method();
	}
	catch (std::exception& error) {
		std::cout << std::endl << error.what() << std::endl;
		return 0;
	}
	return 0;
}
