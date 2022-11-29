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
	std::cout << "\n\n\t   ";
	for (size_t i = 0; i < free.size(); ++i) {
		std::cout << free[i] << "		   ";
	}
	std::cout << "\n";
	for (size_t i = 0; i < basis.size(); ++i) {
		std::cout << basis[i] << "	";
		for (size_t j = 0; j < matr[i].size(); ++j) {
			if (matr[i][j] == -0) {
				matr[i][j] = 0;
			}
			std::cout << std::setw(8) << std::fixed << std::setprecision(3) << round(matr[i][j] * 1000) / 1000 << std::right << "\t";
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
					transformation(matr, index_r_raw, index_r_column);
					std::swap(free[index_r_column], basis[index_r_raw]);
					otr = true;
					break;
				}
			}
			if (otr) {
				printer(matr, basis, free);
				i = 0;
			}
			else {
				throw std::logic_error("No sollution to this task :(");
			}
		}
	}
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
				throw std::logic_error(R"(Infinity number of sollution \(-_-)/)");
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
		if (rezult == "max") {
			std::cout << "W  = " << promezutochniy_result_func << std::endl;
			double g = 0;
			size_t i = 1;
			double sum = 0;
			g = 1 / promezutochniy_result_func;
			std::cout << "\ng  = " << g << std::endl;
			std::cout << "\nOptimal strategy for player 'A' is: \n";
			for (const auto& el : m_result) {
				std::cout << "x" << i << " = " << el.second * g << std::endl;
				sum += el.second * g;
				++i;
			}
			std::cout << "\nChecking: " << sum << " = 1.000";
		}
		if (rezult == "min") {
			std::cout << "Z  = " << promezutochniy_result_func << std::endl;
			double h = 0;
			size_t i = 1;
			double sum = 0;
			h = 1 / promezutochniy_result_func;
			std::cout << "\nh  = " << h << std::endl;
			std::cout << "\nOptimal strategy for player 'B' is: \n";
			for (const auto& el : m_result) {
				std::cout << "y" << i << " = " << el.second * h << std::endl;
				sum += el.second * h;
				++i;
			}
			std::cout << "\nChecking: " << sum << " = 1.000";
		}
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

	std::vector<std::vector<double>> usl = {
	{8, 1, 17, 8, 1},
	{12, 6, 11, 10, 16},
	{4, 19, 11, 15, 2},
	{17, 19, 6, 17, 16}
	};

	std::cout << "Optimal strategy for 'A' player:\n";

	std::vector<std::vector<double>> A_for_a_player(usl[0].size());

	for (size_t i = 0; i < A_for_a_player.size(); ++i) {
		A_for_a_player[i].resize(usl.size());
		for (size_t j = 0; j < A_for_a_player[i].size(); ++j) {
			A_for_a_player[i][j] = (-1) * usl[j][i];
		}
	}

	for (size_t i = 0; i < A_for_a_player.size(); ++i) {
		A_for_a_player[i].insert(A_for_a_player[i].end(), -1);
	}

	std::vector<double> row1(A_for_a_player[0].size());

	for (size_t i = 0; i < row1.size(); ++i) {
		row1[i] = -1;
	}

	A_for_a_player.push_back(row1);

	for (size_t i = 0; i < A_for_a_player.size(); ++i) {
		std::swap(A_for_a_player[i][0], A_for_a_player[i][A_for_a_player.size() - 2]);
		std::swap(A_for_a_player[i][1], A_for_a_player[i][A_for_a_player.size() - 2]);
		std::swap(A_for_a_player[i][2], A_for_a_player[i][A_for_a_player.size() - 2]);
		std::swap(A_for_a_player[i][3], A_for_a_player[i][A_for_a_player.size() - 2]);
	}

	A_for_a_player[A_for_a_player.size() - 1][0] = 0;

	std::vector<std::string> free = { "sv", "u1", "u2", "u3" , "u4" };
	std::vector<std::string> basis = { "u5", "u6", "u7", "u8", "u9", "W" };
	std::string rezult = "max";


	Simplex_tabels ob(basis, free, A_for_a_player, rezult);

	try {
		ob.Simpl_method();
	}
	catch (std::exception& error) {
		std::cout << std::endl << error.what() << std::endl;
		return 0;
	}


	std::cout << "\n\nOptimal strategy for 'B' player:\n";

	std::vector<std::vector<double>>  A_for_b_player(usl.size());

	for (size_t i = 0; i < A_for_b_player.size(); ++i) {
		A_for_b_player[i].resize(usl[i].size());
		for (size_t j = 0; j < A_for_b_player[i].size(); ++j) {
			A_for_b_player[i][j] = usl[i][j];
		}
	}

	for (size_t i = 0; i < A_for_b_player.size(); ++i) {
		A_for_b_player[i].insert(A_for_b_player[i].end(), 1);
	}

	std::vector<double> row(A_for_b_player[0].size());

	for (size_t i = 0; i < row.size(); ++i) {
		row[i] = -1;
	}

	A_for_b_player.push_back(row);

	for (size_t i = 0; i < A_for_b_player.size(); ++i) {
		std::swap(A_for_b_player[i][0], A_for_b_player[i][A_for_b_player.size()]);
		std::swap(A_for_b_player[i][1], A_for_b_player[i][A_for_b_player.size()]);
		std::swap(A_for_b_player[i][2], A_for_b_player[i][A_for_b_player.size()]);
		std::swap(A_for_b_player[i][3], A_for_b_player[i][A_for_b_player.size()]);
		std::swap(A_for_b_player[i][4], A_for_b_player[i][A_for_b_player.size()]);
	}

	A_for_b_player[A_for_b_player.size() - 1][0] = 0;

	std::vector<std::string> freeb = { "sv", "v1", "v2", "v3" , "v4", "v5" };
	std::vector<std::string> basisb = { "v6", "v7", "v8", "v9", "Z" };
	std::string rezultb = "min";

	Simplex_tabels obb(basisb, freeb, A_for_b_player, rezultb);

	try {
		obb.Simpl_method();
	}
	catch (std::exception& error) {
		std::cout << std::endl << error.what() << std::endl;
		return 0;
	}
	return 0;
}
