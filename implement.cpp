#include <bits/stdc++.h>

using namespace std;

// Hàm mục tiêu: ví dụ f(x) = sum(x_i^2)
double fitness(const vector<double>& X) {
    double sum = 0.0;
    for (double xi : X)
        sum += xi * xi;
    return sum;
}

// Giới hạn biến
double clamp(double val, double lb, double ub) {
    if (val < lb) return lb;
    if (val > ub) return ub;
    return val;
}

int main() {
    srand(time(0));

    // Tham số thuật toán
    int dim = 5;             // số chiều
    int N = 20;              // số sói
    int max_iter = 100;      // số vòng lặp tối đa
    double lb = -10.0;       // biên dưới
    double ub = 10.0;        // biên trên

    // Khởi tạo quần thể
    vector<vector<double>> wolves(N, vector<double>(dim));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < dim; j++)
            wolves[i][j] = lb + (ub - lb) * ((double) rand() / RAND_MAX);

    // Khởi tạo alpha, beta, delta
    vector<double> alpha(dim), beta(dim), delta(dim);
    double alpha_score = numeric_limits<double>::infinity();
    double beta_score = numeric_limits<double>::infinity();
    double delta_score = numeric_limits<double>::infinity();

    // Vòng lặp chính
    for (int iter = 0; iter < max_iter; iter++) {
        double a = 2.0 * (1.0 - (double)iter / max_iter); // giảm dần a từ 2 xuống 0

        // Đánh giá quần thể và cập nhật alpha, beta, delta
        for (int i = 0; i < N; i++) {
            double fit = fitness(wolves[i]);
            if (fit < alpha_score) {
                delta_score = beta_score;
                delta = beta;
                beta_score = alpha_score;
                beta = alpha;
                alpha_score = fit;
                alpha = wolves[i];
            }
            else if (fit < beta_score) {
                delta_score = beta_score;
                delta = beta;
                beta_score = fit;
                beta = wolves[i];
            }
            else if (fit < delta_score) {
                delta_score = fit;
                delta = wolves[i];
            }
        }

        // Cập nhật vị trí sói
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < dim; j++) {
                double r1 = (double) rand() / RAND_MAX;
                double r2 = (double) rand() / RAND_MAX;

                double A1 = 2 * a * r1 - a;
                double C1 = 2 * r2;

                double D_alpha = fabs(C1 * alpha[j] - wolves[i][j]);
                double X1 = alpha[j] - A1 * D_alpha;

                r1 = (double) rand() / RAND_MAX;
                r2 = (double) rand() / RAND_MAX;
                double A2 = 2 * a * r1 - a;
                double C2 = 2 * r2;
                double D_beta = fabs(C2 * beta[j] - wolves[i][j]);
                double X2 = beta[j] - A2 * D_beta;

                r1 = (double) rand() / RAND_MAX;
                r2 = (double) rand() / RAND_MAX;
                double A3 = 2 * a * r1 - a;
                double C3 = 2 * r2;
                double D_delta = fabs(C3 * delta[j] - wolves[i][j]);
                double X3 = delta[j] - A3 * D_delta;

                wolves[i][j] = clamp((X1 + X2 + X3) / 3.0, lb, ub);
            }
        }

        // In kết quả mỗi 10 vòng
        if (iter % 10 == 0) {
            cout << "Iter " << iter << ": Best fitness = " << alpha_score << endl;
        }
    }

    // In kết quả cuối cùng
    cout << "\nAlpha solution: ";
    for (double x : alpha) cout << x << " ";
    cout << "\nFitness = " << alpha_score << endl;

    return 0;
}
