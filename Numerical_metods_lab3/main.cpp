#include <boost/numeric/ublas/matrix.hpp>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <json.hpp>
#include "input.hpp"
#include "start_proc.hpp"

#undef max

using matrix = boost::numeric::ublas::matrix<double, boost::numeric::ublas::row_major>;

struct func_t
{
    uint32_t n, m;

    double operator()(int j, int i, bool exact_solution = false)
    {
        assert(j <= n && j >= 0);
        assert(i <= m && i >= 0);

        double x = static_cast<double>(i) / n;
        double y = static_cast<double>(j) / n;

        if (exact_solution)
            return std::pow(x, 3) + std::pow(y, 2) + 3.0;
        return 6 * x + 2;
    }
};

func_t u;

void print_matrix(const matrix& mat, std::ostream& out = std::cout)
{
    for (int i = 0; i < mat.size2(); i++)
    {
        for (int j = 0; j < mat.size1(); j++)
        {
            out << std::right << mat(j, i) << ';';
        }
        out << '\n';
    }

    out << '\n';
}

void init_borders(matrix& v)
{
    for (int i = 0; i < v.size1(); i++)
    {
		v(i, 0) = u(i, 0, true);
        v(i, v.size2() - 1) = u(i, v.size2() - 1, true);
    }

    for (int i = 1; i < v.size2() - 1; i++)
    {
        v(0, i) = u(0, i, true);
		v(v.size1() - 1, i) = u(v.size1() - 1, i, true);
    }
}

matrix solve(const input_t& input)
{
    uint32_t n = input.n;
    uint32_t m = input.m;

    double h = 1.0 / n;
    double k = 1.0 / m;

    double b = n * n; // 1/ h^2
    double c = m * m; // 1/ k^2

    uint64_t max_step = input.max_step;
    double eps = input.eps;

    uint64_t count_step = 0;

	matrix v(n + 1, m + 1);
    
    auto it1 = v.begin1();

	for (; it1 != v.end1(); ++it1)
    {
        std::fill(it1.begin(), it1.end(), 0.0);
    }

    init_borders(v);

    double A = -2.0 * (1 / (h * h) + 1 / (k * k));

    while (count_step++ < max_step)
    {
        double cur_eps = 0;
        for (int i = 1; i < m; i++)
        {
            for (int j = 1; j < n; j++)
            {
				//double tmp_v = -(b * (v(i - 1, j) + v(i + 1, j)) + c * (v(i, j - 1) + v(i, j + 1)) - u(i, j)) / A;
				double tmp_v = -b * (v(j - 1, i) + v(j + 1, i)) - c * (v(j, i - 1) + v(j, i + 1));
                tmp_v += u(j, i);
                tmp_v /= A;

                cur_eps = std::max(std::abs(tmp_v - v(j, i)), cur_eps);

                v(j, i) = tmp_v;
            }
        }

        if (cur_eps < eps)
            break;
    }

    std::cout << "Step count is: " << count_step << '\n';


    return v;
}


matrix get_z(const matrix& actual)
{
    matrix z(actual.size1(), actual.size2());

    for (int i = 0; i < actual.size2(); i++)
    {
        for (int j = 0; j < actual.size1(); j++)
        {
			z(j,i) = u(j, i, true) - actual(j, i);
        }
    }

    return z;
}

matrix get_r(const matrix& v)
{
    matrix r(v.size1(), v.size2());

    uint32_t n = v.size1() - 1;
    uint32_t m = v.size2() - 1;

    auto it1 = r.begin1();

    for (; it1 != r.end1(); ++it1)
    {
        std::fill(it1.begin(), it1.end(), 0.0);
    }

    double b = n * n;
    double c = m * m;
    double A = -2.0 * (b + c);

    for (int i = 1; i < m; i++)
    {
        for (int j = 1; j < n; j++)
        {
            r(j, i) = A * v(j, i) + 
                b * (v(j - 1, i) + v(j + 1, i)) + c * (v(j, i - 1) + v(j, i + 1));

            r(j, i) -= u(j, i);
        }
    }

    return r;
}

double abs(const matrix& mat)
{
    double res = 0;
    for (int i = 0; i < mat.data().size(); i++)
        res = std::max(std::abs(mat.data()[i]), res);

    return res;
}

int main() 
{
    std::cout << std::setprecision(8);

    input_t input;

    std::string input_f_name = "input.json";

    std::ifstream input_f(input_f_name);

    if (!input_f.is_open())
    {
        std::cerr << "Cant open input file\n";
        exit(1);
    }

    nlohmann::json json = nlohmann::json::parse(input_f);

    input_f.close();

    get_input(json, input);

    u.n = input.n;
    u.m = input.m;

    matrix ans = solve(input);

    std::cout << "DONE!, starting output\n";

    std::ofstream out("out.csv");

    out << "v:\n";
    print_matrix(ans, out);

    //double discrepancy = get_z(expected, ans);

    matrix r = get_r(ans);
    matrix z = get_z(ans);

    out << "|r|; " << std::scientific << abs(r) << '\n';
    out << "|z|; " << std::scientific << abs(z) << "\n\n";

    out << std::defaultfloat;

    out << "r:\n";
    print_matrix(r, out);
    out << "z:\n";
    print_matrix(z, out);

    out.close();

    start_proc();
    //std::cout << "Exact solution:\n";
    //print_matrix(expected);
}