//==================================================================//
// Copyright 2017 Clement Ranc
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//==================================================================//
// This file define a class of three dimensional table.
//==================================================================//

#ifndef __TABLE_H_
#define __TABLE_H_

#include<cassert>
#include<iostream>

using namespace std;

class Table {

    // Data
    protected:
    int size1;
    int size2;
    int size3;
    double* tab;

    // Constructors
    public:
    explicit Table(int d1, int d2=1, int d3=1);
    Table(const Table&);
    explicit Table(const char*);
    explicit Table(const char*, const char*, const int);
    explicit Table(const double*, int d1, int d2=1, int d3=1);

    // Destructor
    virtual ~Table();

    // Definitions
    void operator=(const Table&);
    void operator=(double);


    // Reading access to data
    int get_size1() const { return size1; };
    int get_size2() const { return size2; };
    int get_size3() const { return size3; };

    double raw_table(int i) const { return tab[i]; };

    double operator()(int i, int j=0, int k=0) const {
        assert ((i>=0) && (i<size1));
        assert ((j>=0) && (j<size2));
        assert ((k>=0) && (k<size3));
        return tab[(i * size2 + j) * size3 + k];
    };

    // Writing access to data
    double& set(int i, int j=0, int k=0){
        assert ((i>=0) && (i<size1));
        assert ((j>=0) && (j<size2));
        assert ((k>=0) && (k<size3));
        return tab[(i * size2 + j) * size3 + k];
    };

    // Useful functions
    void save(const char* fname) const;

    protected:
    virtual void print(ostream&) const;  // Print the table

    friend ostream& operator<<(ostream&, const Table& );
    friend Table operator-(const Table&);
    friend Table operator+(const Table&, const Table&);
    friend Table operator+(const Table&, double);
    friend Table operator+(double, const Table&);
    friend Table operator-(const Table&, const Table&);
    friend Table operator-(const Table&, double);
    friend Table operator-(double, const Table&);
    friend Table operator*(const Table&, double);
    friend Table operator*(double, const Table&);
    friend Table operator/(const Table&, const Table&);
    friend Table operator/(const Table&, double);
    friend Table operator/(double, const Table&);

    friend Table sqrt(const Table&);
    friend Table pow(const Table&, double);
};

ostream& operator<<(ostream&, const Table& );
Table operator-(const Table&);
Table operator+(const Table&, const Table&);
Table operator+(const Table&, double);
Table operator+(double, const Table&);
Table operator-(const Table&, const Table&);
Table operator-(const Table&, double);
Table operator-(double, const Table&);
Table operator*(const Table&, double);
Table operator*(double, const Table&);
Table operator/(const Table&, const Table&);
Table operator/(const Table&, double);
Table operator/(double, const Table&);

Table sqrt(const Table&);
Table pow(const Table&, double);

#endif
