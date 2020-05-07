#include <iostream>
#include <string>
#include <vector>

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xaxis_slice_iterator.hpp"
#include "xtensor/xadapt.hpp"

#define H5_USE_XTENSOR
#include <highfive/H5Easy.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>


// #include <bits/stdc++.h> 
// #include <boost/algorithm/string.hpp> 

using namespace std;
using namespace xt;

struct Sum
{
   void operator()(const xarray<int>& a) { sum += xt::sum(a)(); }
   int sum {0};
};

int main(int argc, char** argv) {
    xarray<int> a = {{{1, 2, 3, 4},
                  {5, 6, 7, 8},
                  {9, 10, 11, 12}},
                 {{13, 14, 15, 16},
                  {17, 18, 19, 20},
                  {21, 22, 23, 24}}};

    auto iter = axis_slice_begin(a, 0);
    auto end = axis_slice_end(a, 0);

    auto print = [](const xarray<int>& a) {cout << a;};
    std::for_each(iter, end, print);
    
    iter = axis_slice_begin(a, 0);
    Sum s1 = std::for_each(iter, end, Sum());
    cout << endl << "Sum is " << s1.sum << endl;

    const auto& s = a.shape();
    std::copy(s.cbegin(), s.cend(), std::ostream_iterator<double>(std::cout, " "));

    xarray<int> b = {{1, 2, 3, 4},
                    {5, 6, 7, 8}};
    cout << xt::col(b, 0) << endl;
    xt::col(b, 0) = xt::xarray<int>({6, 7});
    cout << xt::col(b, 0) << endl;

    vector<int> tmp {12, 13};
    vector<size_t> shape = {2};
    xt::col(b, 0) = xt::adapt(tmp, shape);
    cout << xt::col(b, 0) << endl;

    return 0;
}