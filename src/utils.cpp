#include "../include/utils.h"

namespace utils
{

    using namespace std;
    namespace ublas = boost::numeric::ublas;

    int GetFullOrderIndex(int r, int s, int order)
    {
        int sum_temp = 0;
        int k;
        for (int sp = 0; sp < s; sp++)
        {
            sum_temp += order + 1 - sp;
        }
        k = sum_temp + r + 1;
        return k - 1; // return the index which is zero based
    }

    std::vector<int> GetVertexIndex(std::vector<int> &element)
    {
        int num_node_in_element = element.size();
        // calculate the order of bases based on the number of nondes in the element
        int num_order = int((sqrt(8 * num_node_in_element + 1) - 3) / 2);
        std::vector<int> vertex_index(3);
        vertex_index[0] = GetFullOrderIndex(0, 0, num_order);
        vertex_index[1] = GetFullOrderIndex(num_order, 0, num_order);
        vertex_index[2] = GetFullOrderIndex(0, num_order, num_order);
        return vertex_index;
    }

    bool SortByColumn0(vector<int> const &v1, vector<int> const &v2)
    {
        return v1[0] < v2[0];
    }

    bool SortByColumn2(vector<int> const &v1, vector<int> const &v2)
    {
        return v1[2] < v2[2];
    }

    template <typename T>
    vector<vector<T>> SliceByRow(vector<vector<T>> const &v, int m, int n)
    {
        auto first = v.cbegin() + m;
        auto last = v.cbegin() + n + 1;
        vector<vector<T>> vec(first, last);
        return vec;
    }
    template vector<vector<int>> SliceByRow<int>(vector<vector<int>> const &v, int m, int n);

    int MissingFrom012(int input1, int input2)
    {
        if ((input1 == 0 && input2 == 1) || (input1 == 1 && input2 == 0))
        {
            return 2;
        }
        else if ((input1 == 0 && input2 == 2) || (input1 == 2 && input2 == 0))
        {
            return 1;
        }
        else if ((input1 == 1 && input2 == 2) || (input1 == 2 && input2 == 1))
        {
            return 0;
        }
        else
        {
            abort();
        }
    }

    template<typename T>
    ublas::vector<T> StdToBoostVector(std::vector<T> std_vec)
    {
        ublas::vector<T> boost_vec (std_vec.size());
        for (int i = 0; i < std_vec.size(); i++)
        {
            boost_vec(i) = std_vec[i];
        }
        return boost_vec;
    }
    template ublas::vector<double> StdToBoostVector<double>(std::vector<double> std_vec);

    double MaxBoostVector(ublas::vector<double> vec)
    {
        double max_value = 0.0;
        for (int i = 0; i < vec.size(); i++)
        {
            if (vec[i] > max_value)
            {
                max_value = vec[i];
            }
        }
        return max_value;
    }

    ublas::matrix<double> Invert22Matrix(ublas::matrix<double> mat_input)
    {
        double a, b, c, d;
        ublas::matrix<double> mat_inversed(2, 2, 0.0);
        a = mat_input(0, 0); b = mat_input(0, 1);
        c = mat_input(1, 0); d = mat_input(1, 1);
        double det = a * d - b * c;
        mat_inversed(0, 0) = d / det;
        mat_inversed(0, 1) = -b / det;
        mat_inversed(1, 0) = -c / det;
        mat_inversed(1, 1) = a /det;
        return mat_inversed;
    }
} // namespace utils
