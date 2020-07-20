//Source: https://gist.githubusercontent.com/hallazzang/4e6abbb05ff2d3e168a87cf10691c4fb/raw/9c1461149b632777d830734e34264c5e96f03004/Matrix.java
//Modified by Sanket S. Desai / 04-08-2018
package jabix;

public class IntegerMatrix {
    int[][] data = null;
    int rows = 0, cols = 0;

    public IntegerMatrix(){
      int[][] data=new int[0][0];
    }
    public IntegerMatrix(int rows, int cols) {
        this.data = new int[rows][cols];
        for(int i=0; i < rows; i++){
          for(int j=0; j < cols; j++){
            this.setElement(i,j,0);
          }
        }
        this.rows = rows;
        this.cols = cols;
    }

    public IntegerMatrix(int[][] data) {
        this.data = data.clone();
        rows = this.data.length;
        cols = this.data[0].length;
    }
    void initialize(int rws, int cls){
      this.data=new int[rws][cls];
      for(int i=0; i < rws; i++){
        for(int j=0; j < cls; j++){
          this.setElement(i,j,0);
        }
      }
    }
    int getRowSum(int rw){
      int rsum=0;
      for(int i=0; i < this.cols; i++ ){
        rsum+=this.data[rw][i];
      }
      return rsum;
    }
    void setNumberOfRows(int r)
    {
      this.rows=r;
    }
    void setNumberOfColumns(int c){
      this.cols=c;
    }
    int getColumnSum(int cl){
      int csum=0;
      for(int i=0; i< this.rows; i++){
        csum+=this.data[i][cl];
      }
      return csum;
    }
    boolean isSquare() {
        return rows == cols;
    }

    int getElement(int i, int j){
      return this.data[i][j];
    }

    void setElement(int i, int j, int value){
      this.data[i][j]=value;
    }
    int getNumberOfRows(){
      return this.rows;
    }
    int getNumberOfColumns(){
      return this.cols;
    }
    void display() {
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

    IntegerMatrix transpose() {
        IntegerMatrix result = new IntegerMatrix(cols, rows);

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                result.data[col][row] = data[row][col];
            }
        }

        return result;
    }

    // Note: exclude_row and exclude_col starts from 1
    static IntegerMatrix subMatrix(IntegerMatrix matrix, int exclude_row, int exclude_col) {
        IntegerMatrix result = new IntegerMatrix(matrix.rows - 1, matrix.cols - 1);

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

    int determinant() {
        if (rows != cols) {
            //return Integer.NaN;
            return -999999999;
        }
        else {
            return _determinant(this);
        }
    }

    private int _determinant(IntegerMatrix matrix) {
        if (matrix.cols == 1) {
            return matrix.data[0][0];
        }
        else if (matrix.cols == 2) {
            return (matrix.data[0][0] * matrix.data[1][1] -
                    matrix.data[0][1] * matrix.data[1][0]);
        }
        else {
            int result = 0;

            for (int col = 0; col < matrix.cols; ++col) {
                IntegerMatrix sub = subMatrix(matrix, 1, col + 1);

                result += (Math.pow(-1, 1 + col + 1) *
                           matrix.data[0][col] * _determinant(sub));
            }

            return result;
        }
    }

    IntegerMatrix inverse() {
        int det = determinant();

        if (rows != cols || det == 0.0) {
            return null;
        }
        else {
            IntegerMatrix result = new IntegerMatrix(rows, cols);
            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    IntegerMatrix sub = subMatrix(this, row + 1, col + 1);

                    result.data[col][row] = (1 / det *
                                            (int)( Math.pow(-1, row + col)) *
                                             _determinant(sub));
                }
            }

            return result;
        }
    }
}
