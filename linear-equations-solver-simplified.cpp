#include <iostream>

using namespace std;

/*
 * Solve the linear system A x = b by Gauss-Jordan elimination, where A is an
 * rows-by-cols coefficient matrix and b is the constant vector.
 *
 * Input
 *     rows cols
 *     then rows lines, each holding cols coefficients followed by the constant
 *     of that equation
 *
 * Output, one of the three possibilities
 *     "The solution exists and is unique." + "Solution: ( ... )T"
 *         A is square and nonsingular, so the reduced row echelon form itself
 *         is "x = the constant column"
 *     "There are infinitely many solutions." + "Solution: ( p )T + xF1( w1 )T + ..."
 *         p is one particular solution and every wi is a vector of the null
 *         space, so every solution is p plus a linear combination of the wi
 *     "There is no solution."
 *         The system is inconsistent
 *
 * Complexity
 *     O(rows^2 * cols) arithmetic operations
 */

const int MAXSIZE = 102;      // largest accepted number of rows and of columns

const double TOLERANCE = 1e-7; // an entry is treated as zero below this value

const int FEASIBLE = 0;    // exit code: the system is consistent
const int INFEASIBLE = 1;  // exit code: the system is inconsistent

double A[MAXSIZE][MAXSIZE]; // coefficient matrix
double b[MAXSIZE];          // constant vector
int pivotColumn[MAXSIZE];   // column of the leading entry of every nonzero row

/*
 * Subtract a multiple of the pivot row from a row below (or above) it, so that
 * the entry of that row in the pivot column becomes exactly zero.
 *
 * The elementary row operation used here is
 *     row[target] <- row[target] - row[target][pivot] * row[pivot]
 * whose right hand side follows the same rule.
 */
void eliminate(int pivotRow, int targetRow, int cols)
{
    for (int col = pivotColumn[pivotRow] + 1; col <= cols; col++)
    {
        A[targetRow][col] -= A[pivotRow][col] * A[targetRow][pivotColumn[pivotRow]];
    }
    b[targetRow] -= b[pivotRow] * A[targetRow][pivotColumn[pivotRow]];
    A[targetRow][pivotColumn[pivotRow]] = 0;
}

int main()
{
    // Define and input
    int rows, cols; // number of rows and columns of the coefficient matrix
    cin >> rows >> cols;
    for (int row = 1; row <= rows; row++)
    {
        for (int col = 1; col <= cols; col++)
        {
            cin >> A[row][col];
        }
        cin >> b[row];
    }

    // Forward elimination
    // Process one pivot after another. At the end of every round the pivot of
    // the current row is 1 and every entry below it is 0, so the matrix is
    // turned into row echelon form.
    int rank = 0; // number of nonzero rows, which is also the rank of A
    for (int row = 1; row <= rows && rank < cols; row++)
    {
        int pivot = row; // column of the leading entry of this round

        // Find the first column, from the current one on, that still has a
        // nonzero entry below, and take the entry of biggest absolute value as
        // the pivot to limit the growth of the rounding error
        int bestRow = row; // row that holds that biggest entry
        bool found = false; // existence of such a column
        for (pivot = row; pivot <= cols; pivot++)
        {
            for (int candidate = row; candidate <= rows; candidate++)
            {
                if (abs(A[candidate][pivot]) > abs(A[bestRow][pivot]))
                {
                    bestRow = candidate;
                }
            }
            if (abs(A[bestRow][pivot]) > TOLERANCE)
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            // No pivot is left, so every remaining row is a zero row
            break;
        }
        pivotColumn[++rank] = pivot; // record the pivot column of this row

        // Deal with the row with the biggest first element
        if (row != bestRow)
        {
            swap(A[bestRow], A[row]);
            swap(b[bestRow], b[row]);
        }
        for (int col = pivot + 1; col <= cols; col++)
        {
            A[row][col] /= A[row][pivot];
        }
        b[row] /= A[row][pivot];
        A[row][pivot] = 1;

        // Eliminate other rows
        for (int target = row + 1; target <= rows; target++)
        {
            eliminate(row, target, cols);
        }
    }

    // Judge the solution
    // Every row below the nonzero rows is a zero row, so such a row means an
    // equation "0 = b[row]"; the system is inconsistent iff any of those
    // constants is nonzero
    for (int row = rank + 1; row <= rows; row++)
    {
        if (abs(b[row]) > TOLERANCE)
        {
            cout << "There is no solution.";
            return INFEASIBLE;
        }
    }

    // Back subtitution
    // Clear the entries above every pivot as well, which reduces the matrix
    // further into reduced row echelon form. The forward elimination has
    // already made every pivot 1, so the pivot column of a row is known and
    // there is nothing to search for here.
    for (int row = rank; row >= 1; row--)
    {
        // Eliminate other rows
        for (int target = 1; target <= row - 1; target++)
        {
            eliminate(row, target, cols);
        }
    }

    if (rank == cols)
    {
        // Situation 2: unique solution
        // Every column is a pivot column, so the reduced system reads x = b
        cout << "The solution exists and is unique." << endl;
        cout << "Solution: ( ";
        for (int row = 1; row <= cols; row++)
        {
            cout << b[row] << " ";
        }
        cout << ")T";
    }
    else
    {
        // Situation 3: infinitely many solutions
        cout << "There are infinitely many solutions." << endl;

        // A column is a pivot column iff some row has its first nonzero entry in it
        bool isPivot[MAXSIZE];
        for (int col = 1; col <= cols; col++)
        {
            isPivot[col] = false;
        }
        for (int row = 1; row <= rank; row++)
        {
            isPivot[pivotColumn[row]] = true; // the pivot columns are already known
        }

        // Find free columns (columns with free variables)
        // A free column carries no pivot, so its variable can be chosen freely
        int freeColumn[MAXSIZE];   // free columns
        bool isFree[MAXSIZE];      // judge each column is or is not free
        int freeCount = 0;         // number of free columns
        for (int col = 1; col <= cols; col++)
        {
            if (!isPivot[col])
            {
                isFree[col] = true;
                freeCount++;
                freeColumn[freeCount] = col;
            }
            else
            {
                isFree[col] = false;
            }
        }

        // Output
        cout << "Solution: ( ";

        // Output constant part;
        // The pivot variables run through the rows in order, while the free
        // variables are set to zero, which leaves the constant column as the
        // particular solution
        int pivotRow = 0; // row that holds the pivot of the current column
        for (int col = 1; col <= cols; col++)
        {
            if (isFree[col])
            {
                cout << "0 ";
            }
            else
            {
                pivotRow++;
                cout << b[pivotRow] << " ";
            }
        }
        cout << ")T ";

        // Output free part;
        // Setting one free variable to 1 and the others to 0 makes that column
        // of the reduced matrix the negative of a null space basis vector
        for (int freeIndex = 1; freeIndex <= freeCount; freeIndex++)
        {
            pivotRow = 0;
            cout << "+ x" << freeColumn[freeIndex] << " ( ";
            for (int col = 1; col <= cols; col++)
            {
                if (isFree[col])
                {
                    if (col == freeColumn[freeIndex])
                    {
                        cout << "1 ";
                    }
                    else
                    {
                        cout << "0 ";
                    }
                }
                else
                {
                    pivotRow++;
                    cout << -A[pivotRow][freeColumn[freeIndex]] << " ";
                }
            }
            cout << ")T ";
        }
    }
    // cin.get();
    // cin.get();
    return FEASIBLE;
}
