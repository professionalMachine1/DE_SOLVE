#include "SDE.h"
#include <sciplot/sciplot.hpp>

namespace sp = sciplot;

void function(double t_, const gsl_vector* x_, gsl_vector* res, const parametrs& p)
{
    res->data[0] = p.a - sin(x_->data[0]);
}

int main()
{
    parametrs p = { 1.5 };
    std::vector<std::vector<double>> res(1 + 1);
    DE_Solve solve(function, pow(10, -6), pow(10, -3));
    // ---
    gsl_vector* x0 = gsl_vector_alloc(1);
    gsl_vector_set(x0, 0, M_PI);

    // Считаем
    solve.set_parametrs(p);
    solve.solve(0, 10 * M_PI, x0, res);

    // Строим
    sp::Plot2D plot;
    plot.drawCurve(res[0], res[1]).label("DE");
    sp::Figure figure = { {plot} };
    sp::Canvas canvas = { {figure} };
    canvas.show();
}
