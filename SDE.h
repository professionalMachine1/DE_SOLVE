#pragma once

#ifndef SDE_H
#define SDE_H

#include <list>
#include <vector>
#include <iostream>
#include <gsl/gsl_linalg.h>

void gsl_vector_abs(gsl_vector* v)
{
    for (size_t i = 0; i < v->size; ++i) {
        v->data[i] = fabs(v->data[i]);
    }
}

void gsl_vector_sub_new(const gsl_vector* a, const gsl_vector* b, gsl_vector* result)
{
    gsl_vector_memcpy(result, a);
    gsl_vector_sub(result, b);
}

void push_vector(const double a_, const gsl_vector* v, std::vector<std::vector<double>>& res)
{
    res[0].push_back(a_);
    for (size_t i = 1; i < res.size(); ++i) {
        res[i].push_back(gsl_vector_get(v, i - 1));
    }
}

void pop_vector(std::vector<std::vector<double>>& res)
{
    for (size_t i = 0; i < res.size(); ++i) {
        res[i].pop_back();
    }
}

struct parametrs
{
    double a = 1;
};

class DE_Solve
{
private:
    parametrs p;
    double step, eps, t;
    gsl_vector* x, * spec_x;
    void(*function)(double t_, const gsl_vector* x_, gsl_vector* res, const parametrs& p_);
public:
    DE_Solve();
    DE_Solve(void(*function_)(double, const gsl_vector*, gsl_vector*, const parametrs&),
        double eps_ = pow(10, -6), double step_ = pow(10, -3));
    void set_parametrs(const parametrs& p_) { p = p_; }
    void euler_method(double t_, double step_, gsl_vector* x_);
    void solve_const_step(double t_min, double t_max, gsl_vector* x0,
        std::vector<std::vector<double>>& vals);
    void solve(double t_min, double t_max, gsl_vector* x0, std::vector<std::vector<double>>& vals);

    ~DE_Solve();
};

DE_Solve::DE_Solve()
{
    step = 0, eps = 0, t = 0;
    x = gsl_vector_alloc(0);
    spec_x = gsl_vector_alloc(0);
    function = NULL;
}

DE_Solve::DE_Solve(void(*function_)(double, const gsl_vector*, gsl_vector*, const parametrs&), double eps_, double step_ )
{
    step = step_, eps = eps_, t = 0;
    x = gsl_vector_alloc(0);
    spec_x = gsl_vector_alloc(0);
    function = function_;
}

void DE_Solve::euler_method(double t_, double step_, gsl_vector* x_)
{
    function(t_, x_, spec_x, p); // spec_x = function(t_, x_)
    gsl_blas_daxpy(step_, spec_x, x_); // x_ = x_ + step_ * spec_x
}

void DE_Solve::solve_const_step(double t_min, double t_max, gsl_vector* x0, 
    std::vector<std::vector<double>>& vals)
{
    double tprev = t_min;
    size_t size = x0->size;
    std::vector<double> cur_val(size + 1);
    gsl_vector* xprev = gsl_vector_alloc(size);

    gsl_vector_free(x);
    gsl_vector_free(spec_x);
    // ---
    x = gsl_vector_alloc(size);
    spec_x = gsl_vector_alloc(size);
    // ---
    gsl_vector_memcpy(x, x0);
    gsl_vector_memcpy(xprev, x0);

    t = t_min;
    while (t < t_max)
    {
        tprev = t;
        gsl_vector_memcpy(xprev, x);

        euler_method(t, step, x);
        t += step;

        push_vector(t, x, vals);
    }
    pop_vector(vals);
    euler_method(tprev, t_max - tprev, xprev);
    // ---
    push_vector(t_max, xprev, vals);
}

void DE_Solve::solve(double t_min, double t_max, gsl_vector* x0, std::vector<std::vector<double>>& vals)
{
    size_t size = x0->size, p = 1; // p - method order
    std::vector<double> cur_val(size + 1);
    double max_S, diff_corr = 1 / (pow(2, p) - 1), spec_eps = eps / pow(2, p + 1), tprev = t_min;
    gsl_vector* S, * X, * xprev;

    gsl_vector_free(x);
    gsl_vector_free(spec_x);
    // ---
    x = gsl_vector_alloc(size);
    X = gsl_vector_alloc(size);
    S = gsl_vector_alloc(size);
    xprev = gsl_vector_alloc(size);
    spec_x = gsl_vector_alloc(size);
    // ---
    gsl_vector_memcpy(x, x0);
    gsl_vector_memcpy(xprev, x0);

    t = t_min;
    while (t < t_max)
    {
        tprev = t;
        gsl_vector_memcpy(X, x); // X = x
        gsl_vector_memcpy(xprev, x); // xprev = x

        // Âûחûגאול לועמה Ýיכונא
        euler_method(t, step, x);
        euler_method(t, step * 0.5, X);
        euler_method(t + step * 0.5, step * 0.5, X);

        gsl_vector_sub_new(x, X, S); // S = x - X
        gsl_vector_abs(S);
        max_S = gsl_vector_max(S);

        // ---
        if (max_S > eps)
        {
            gsl_vector_memcpy(x, xprev);
            step *= 0.5;
        }
        else
        {
            t += step;
            // ---
            push_vector(t, x, vals);

            if (max_S < spec_eps)
                step *= 2;
        }
    }
    pop_vector(vals);
    euler_method(tprev, t_max - tprev, xprev);
    // ---
    push_vector(t_max, xprev, vals);

    gsl_vector_free(S);
    gsl_vector_free(X);
    gsl_vector_free(xprev);
}

DE_Solve::~DE_Solve()
{
    gsl_vector_free(x);
    gsl_vector_free(spec_x);
}

#endif
