// Validation adapters call the repository's unmodified numerical C++ sources.
#include "mex.h"
#include "svm.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

void segment_gateway(int, mxArray **, int, const mxArray **);
void circum_gateway(int, mxArray **, int, const mxArray **);

extern "C" int segment(const double *p, int n, double width)
{
    mxArray a(n, 2), w(1, 1);
    std::copy(p, p + 2 * n, a.values.begin());
    w.values[0] = width;
    const mxArray *inputs[] = {&a, &w};
    mxArray *outputs[2] = {};
    segment_gateway(2, outputs, 2, inputs);
    int end = static_cast<int>(outputs[1]->values[0]);
    delete outputs[0];
    delete outputs[1];
    return end;
}

extern "C" void circum(const double *p, double *out)
{
    mxArray a(1, 3), b(1, 3), c(1, 3);
    std::copy(p, p + 3, a.values.begin());
    std::copy(p + 3, p + 6, b.values.begin());
    std::copy(p + 6, p + 9, c.values.begin());
    const mxArray *inputs[] = {&a, &b, &c};
    // Numerical bridge also exercises the production single-output call.
    mxArray *outputs[1] = {};
    circum_gateway(1, outputs, 3, inputs);
    std::copy(outputs[0]->values.begin(), outputs[0]->values.end(), out);
    for (auto *v : outputs)
        delete v;
}

static void quiet(const char *) {}
extern "C" int classify(const double *x, const double *y, int n, int d, const double *test, int nt,
                        double *pred, double *probs)
{
    std::vector<std::vector<svm_node>> nodes(n, std::vector<svm_node>(d + 1));
    std::vector<svm_node *> ptrs(n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < d; ++j)
            nodes[i][j] = {j + 1, x[i * d + j]};
        nodes[i][d] = {-1, 0};
        ptrs[i] = nodes[i].data();
    }
    svm_problem problem{n, const_cast<double *>(y), ptrs.data()};
    svm_parameter param{};
    param.svm_type = C_SVC;
    param.kernel_type = LINEAR;
    param.degree = 3;
    param.gamma = 2.79e-4;
    param.cache_size = 100;
    param.eps = 0.001;
    param.C = 1;
    param.nu = 0.5;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 1;
    if (svm_check_parameter(&problem, &param))
        return -1;
    svm_set_print_string_function(quiet);
    std::srand(240925);
    auto *model = svm_train(&problem, &param);
    int classes = svm_get_nr_class(model);
    std::vector<double> probabilities(classes);
    std::vector<svm_node> row(d + 1);
    for (int i = 0; i < nt; ++i)
    {
        for (int j = 0; j < d; ++j)
            row[j] = {j + 1, test[i * d + j]};
        row[d] = {-1, 0};
        pred[i] = svm_predict_probability(model, row.data(), probabilities.data());
        probs[i] = 0;
        for (double value : probabilities)
            probs[i] += value;
    }
    svm_free_and_destroy_model(&model);
    return classes;
}

#ifdef VALIDATE_MAIN
template <int N> int check_outputs()
{
    mxArray a(1, 3), b(1, 3), c(1, 3);
    a.values = {3, 4, 5};
    b.values = {5, 4, 5};
    c.values = {3, 6, 5};
    const mxArray *in[] = {&a, &b, &c};
    mxArray *out[N] = {}; // Exact caller-sized storage under ASan.
    circum_gateway(N, out, 3, in);
    bool ok = out[0]->values == std::vector<double>({1, 1, 0});
    if (N >= 2) ok = ok && out[1]->values[0] == 0.5;
    if (N >= 3) ok = ok && out[2]->values[0] == 0.5;
    for (auto *value : out) delete value;
    return ok ? 0 : 4;
}
int main(int argc, char **argv)
{
    if (argc > 1 && std::strcmp(argv[1], "output-contract") == 0)
        return check_outputs<1>() || check_outputs<2>() || check_outputs<3>();
    double tri[] = {0, 0, 0, 1, 0, 0, 0, 1, 0}, center[3];
    circum(tri, center);
    if (std::abs(center[0] - 0.5) > 1e-12 || std::abs(center[1] - 0.5) > 1e-12)
        return 2;
    double points[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
    return segment(points, 5, 5) == 5 ? 0 : 3;
}
#endif
