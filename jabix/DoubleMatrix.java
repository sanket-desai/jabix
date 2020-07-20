//Source: https://gist.githubusercontent.com/hallazzang/4e6abbb05ff2d3e168a87cf10691c4fb/raw/9c1461149b632777d830734e34264c5e96f03004/Matrix.java
package jabix;
public class DoubleMatrix {
    private double[][] data = null;
    private int rows = 0, cols = 0;

    public DoubleMatrix(int rows, int cols) {
        data = new double[rows][cols];
        this.rows = rows;
        this.cols = cols;
    }

    public DoubleMatrix(double[][] data) {
        this.data = data.clone();
        rows = this.data.length;
        cols = this.data[0].length;
    }

    public boolean isSquare() {
        return rows == cols;
    }

    public void display() {
        System.out.print("[");
        for (int row = 0; row < rows; ++row) {
            if (row != 0) {
                System.out.print(" ");
            }

            System.out.print("[");

            for (int col = 0; col < cols; ++col) {
                System.out.printf("%8.3f", data[row][col]);

                if (col != cols - 1) {
                    System.out.print(" ");
                }
            }

            System.out.print("]");

            if (row == rows - 1) {
                System.out.print("]");
            }

            System.out.println();
        }
    }

    public DoubleMatrix transpose() {
        DoubleMatrix result = new DoubleMatrix(cols, rows);

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                result.data[col][row] = data[row][col];
            }
        }

        return result;
    }

    // Note: exclude_row and exclude_col starts from 1
    public static DoubleMatrix subMatrix(DoubleMatrix matrix, int exclude_row, int exclude_col) {
        DoubleMatrix result = new DoubleMatrix(matrix.rows - 1, matrix.cols - 1);

        for (int row = 0, p = 0; row < matrix.rows; ++row) {
            if (row != exclude_row - 1) {
                for (int col = 0, q = 0; col < matrix.cols; ++col) {
                    if (col != exclude_col - 1) {
                        result.data[p][q] = matrix.data[row][col];

                        ++q;
                    }
                }

                ++p;
            }
        }

        return result;
    }

    public double determinant() {
        if (rows != cols) {
            return Double.NaN;
        }
        else {
            return _determinant(this);
        }
    }

    private double _determinant(DoubleMatrix matrix) {
        if (matrix.cols == 1) {
            return matrix.data[0][0];
        }
        else if (matrix.cols == 2) {
            return (matrix.data[0][0] * matrix.data[1][1] -
                    matrix.data[0][1] * matrix.data[1][0]);
        }
        else {
            double result = 0.0;

            for (int col = 0; col < matrix.cols; ++col) {
                DoubleMatrix sub = subMatrix(matrix, 1, col + 1);

                result += (Math.pow(-1, 1 + col + 1) *
                           matrix.data[0][col] * _determinant(sub));
            }

            return result;
        }
    }

    public DoubleMatrix inverse() {
        double det = determinant();

        if (rows != cols || det == 0.0) {
            return null;
        }
        else {
            DoubleMatrix result = new DoubleMatrix(rows, cols);

            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    DoubleMatrix sub = subMatrix(this, row + 1, col + 1);

                    result.data[col][row] = (1.0 / det *
                                             Math.pow(-1, row + col) *
                                             _determinant(sub));
                }
            }

            return result;
        }
    }
}
