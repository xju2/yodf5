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
#include <map>


// #include <bits/stdc++.h> 
// #include <boost/algorithm/string.hpp> 

using namespace std;
using namespace xt;

struct Sum
{
   void operator()(const xarray<int>& a) { sum += xt::sum(a)(); }
   int sum {0};
};

template <typename T> 
void print_shape(xt::xarray<T>& a){
    const auto& s = a.shape();
    std::copy(s.cbegin(), s.cend(), std::ostream_iterator<double>(std::cout, " "));
}

int main(int argc, char** argv) {
    xarray<int> a = {{{1, 2, 3, 4},
                  {5, 6, 7, 8},
                  {9, 10, 11, 12}},
                 {{13, 14, 15, 16},
                  {17, 18, 19, 20},
                  {21, 22, 23, 24}}};

    auto iter = axis_slice_begin(a, 0);
    auto end = axis_slice_end(a, 0);

    cout << "print 1 " << endl;
    auto print = [](const xarray<int>& a) {cout << a;};
    std::for_each(iter, end, print);
    
    // iter = axis_slice_begin(a, 0);
    // Sum s1 = std::for_each(iter, end, Sum());
    // cout << endl << "Sum is " << s1.sum << endl;

    cout << "print 2" << endl;
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

	cout << "begin experiment 3" << endl;
	std::vector<size_t> shape3d = {2, 4, 4};
	xt::xarray<int> array3d(shape3d, 0);
    // auto v2 = xt::view(array3d, xt::all() , xt::all(), 1);
    // auto v2 = xt::view(array3d, 1, xt::all(), xt::range(1, 3));

	// auto v2 = xt::view(array3d, xt::all(), xt::all(), 1); // failed
	auto v2 = xt::view(array3d, 0, xt::all(), xt::all()); // OK
    // auto v2 = xt::view(array3d, xt::all(), xt::all(), xt::range(0, 1)); // OK
    // auto v2 = xt::view(array3d, xt::all(), xt::range(1, 3), xt::all()); // OK
    // auto v2 = xt::view(array3d, 0, xt::all(), xt::range(1, 3));
    cout << v2.dimension() << endl;
	v2(0, 0) = 1;
	v2 = 3;
    cout << v2 << endl;
	std::vector<size_t> shape2d = {4, 2};
	xt::xarray<int> array2d(shape2d, 0);
	array2d = v2;
    cout << array2d << endl;
	cout << "end experiment 3" << endl;



    xarray<double> a1 = {{0., 1., 2.}, {3., 4., 5.}};
    double b1 = 1.2;
    auto tr = view(a1, 0, all());
    tr = b1;
    a1(0, 0) = 12.0;
    cout << tr << endl;

	// a map of xarray
    using arr2d_t = xt::xarray<double>; 
    std::map<std::string, arr2d_t > data_map;
    vector<string> field_keys = {"Sumw2", "Sumw"};
    long unsigned int n_bins = 2;
    long unsigned int n_variations = 3;
    for(auto field_key: field_keys) {
        data_map[field_key] = xt::xarray<double>::from_shape({n_bins, n_variations});
    }
	cout << data_map["Sumw2"] << endl;
	data_map["Sumw2"](0, 0) = 12.0;
	cout << data_map["Sumw2"] << endl;

	cout << "start 5 " << endl;
	int icc = 12;
	cout << icc++ << endl;
	cout << icc << endl;
	cout << ++icc << endl;
	cout << icc << endl;
	cout << "end 5 " << endl;

    return 0;
}
