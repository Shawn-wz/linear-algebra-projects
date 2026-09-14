#include <iomanip>
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
 *     "Invalid input."
 *         The first line does not declare a usable size, or the coefficients
 *         are missing or are not numbers
 *
 * Complexity
 *     O(rows^2 * cols) arithmetic operations
 *
 * Numerical notes
 *     The zero test is relative to the magnitude of the data instead of being
 *     absolute, so multiplying the whole system by a positive constant does
 *     not change the answer. Partial pivoting keeps the rounding error of the
 *     elimination in check.
 */

const int MAXSIZE = 102; // largest accepted number of rows and of columns

const double RELATIVE_TOLERANCE = 1e-10; // a value is zero below this times the scale

const int FEASIBLE = 0;   // exit code: the system is consistent
const int INFEASIBLE = 1; // exit code: the system is inconsistent
const int BAD_INPUT = 2;  // exit code: the input could not be read

const int DIGITS = 15; // significant digits printed for every number

/*
 * The whole system, so that the helper functions below can work on it without
 * reaching for variables that live in main.
 */
struct System
{
    int rows = 0;                  // number of equations
    int cols = 0;                  // number of unknowns
    double A[MAXSIZE][MAXSIZE]{};  // coefficient matrix
    double b[MAXSIZE]{};           // constant vector
};

/*
 * Read a whole number. Returns false once the stream is no longer readable,
 * which happens for a missing value or for something that is not a number.
 */
bool readInt(int &value)
{
    cin >> value;
    return !cin.fail();
}

/*
 * Read a real number. Returns false once the stream is no longer readable.
 */
bool readDouble(double &value)
{
    cin >> value;
    return !cin.fail();
}

/*
 * Read the whole system. Returns false, leaving nothing half read, as soon as
 * a value is missing or is not a number.
 */
bool readSystem(System &s)
{
    if (!readInt(s.rows) || !readInt(s.cols))
    {
        return false;
    }
    if (s.rows < 1 || s.rows > MAXSIZE || s.cols < 1 || s.cols > MAXSIZE)
    {
        return false;
    }
    for (int row = 1; row <= s.rows; row++)
    {
        for (int col = 1; col <= s.cols; col++)
        {
            if (!readDouble(s.A[row][col]))
            {
                return false;
            }
        }
        if (!readDouble(s.b[row]))
        {
            return false;
        }
    }
    return true;
}

/*
 * The magnitude of the largest entry of the matrix and of the constant vector.
 * An entry much smaller than this is rounding noise rather than data.
 */
double matrixScale(const System &s)
{
    double scale = 0.0;
    for (int row = 1; row <= s.rows; row++)
    {
        for (int col = 1; col <= s.cols; col++)
        {
            if (abs(s.A[row][col]) > scale)
            {
                scale = abs(s.A[row][col]);
            }
        }
        if (abs(s.b[row]) > scale)
        {
            scale = abs(s.b[row]);
        }
    }
    return scale;
}

/*
 * The zero test used everywhere below. Deriving it from the scale of the data
 * instead of fixing it makes the whole program independent of the unit the
 * coefficients are measured in: a system and the same system multiplied by any
 * positive constant get exactly the same treatment.
 *
 * There is deliberately no absolute lower bound here. A floor would make the
 * program treat genuinely small data, such as coefficients of 1e-11, as zero.
 */
double zeroTolerance(const System &s)
{
    return RELATIVE_TOLERANCE * matrixScale(s);
}

/*
 * Subtract a multiple of the pivot row from another row, so that the entry of
 * that row in the pivot column becomes exactly zero.
 *
 * The elementary row operation used here is
 *     row[target] <- row[target] - row[target][column] * row[pivotRow]
 * and the constant vector follows the same rule.
 */
void eliminate(System &s, int pivotRow, int targetRow, int column)
{
    for (int col = column + 1; col <= s.cols; col++)
    {
        s.A[targetRow][col] -= s.A[pivotRow][col] * s.A[targetRow][column];
    }
    s.b[targetRow] -= s.b[pivotRow] * s.A[targetRow][column];
    s.A[targetRow][column] = 0;
}

int main()
{
    // Define and input
    System s;
    if (!readSystem(s))
    {
        cout << "Invalid input.";
        return BAD_INPUT;
    }

    // Every number below is printed with enough digits to be read back
    cout << setprecision(DIGITS);

    double tolerance = zeroTolerance(s);
    int pivotOfRow[MAXSIZE]; // column of the leading entry of every nonzero row

    // Forward elimination
    // Process one pivot after another. At the end of every round the pivot of
    // the current row is 1 and every entry below it is 0, so the matrix is
    // turned into row echelon form.
    int rank = 0; // number of nonzero rows, which is also the rank of A
    for (int row = 1; row <= s.rows && rank < s.cols; row++)
    {
        int pivot = row; // column of the leading entry of this round

        // Find the first column, from the current one on, that still has a
        // nonzero entry below, and take the entry of biggest absolute value as
        // the pivot to limit the growth of the rounding error
        int bestRow = row;  // row that holds that biggest entry
        bool found = false; // existence of such a column
        for (pivot = row; pivot <= s.cols; pivot++)
        {
            for (int candidate = row; candidate <= s.rows; candidate++)
            {
                if (abs(s.A[candidate][pivot]) > abs(s.A[bestRow][pivot]))
                {
                    bestRow = candidate;
                }
            }
            if (abs(s.A[bestRow][pivot]) > tolerance)
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
        pivotOfRow[++rank] = pivot; // record the pivot column of this row

        // Deal with the row with the biggest first element
        if (row != bestRow)
        {
            swap(s.A[bestRow], s.A[row]);
            swap(s.b[bestRow], s.b[row]);
        }
        for (int col = pivot + 1; col <= s.cols; col++)
        {
            s.A[row][col] /= s.A[row][pivot];
        }
        s.b[row] /= s.A[row][pivot];
        s.A[row][pivot] = 1;

        // Eliminate other rows
        for (int target = row + 1; target <= s.rows; target++)
        {
            eliminate(s, row, target, pivot);
        }
    }

    // Judge the solution
    // Every row below the nonzero rows is a zero row, so such a row means an
    // equation "0 = b[row]"; the system is inconsistent iff any of those
    // constants is nonzero
    for (int row = rank + 1; row <= s.rows; row++)
    {
        if (abs(s.b[row]) > tolerance)
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
            eliminate(s, row, target, pivotOfRow[row]);
        }
    }

    if (rank == s.cols)
    {
        // Situation 2: unique solution
        // Every column is a pivot column, so the reduced system reads x = b
        cout << "The solution exists and is unique." << endl;
        cout << "Solution: ( ";
        for (int row = 1; row <= s.cols; row++)
        {
            cout << s.b[row] << " ";
        }
        cout << ")T";
    }
    else
    {
        // Situation 3: infinitely many solutions
        cout << "There are infinitely many solutions." << endl;

        // A column is a pivot column iff some row has its first nonzero entry in it
        bool isPivot[MAXSIZE];
        for (int col = 1; col <= s.cols; col++)
        {
            isPivot[col] = false;
        }
        for (int row = 1; row <= rank; row++)
        {
            isPivot[pivotOfRow[row]] = true; // the pivot columns are already known
        }

        // Find free columns (columns with free variables)
        // A free column carries no pivot, so its variable can be chosen freely
        int freeColumn[MAXSIZE]; // free columns
        bool isFree[MAXSIZE];    // judge each column is or is not free
        int freeCount = 0;       // number of free columns
        for (int col = 1; col <= s.cols; col++)
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
        for (int col = 1; col <= s.cols; col++)
        {
            if (isFree[col])
            {
                cout << "0 ";
            }
            else
            {
                pivotRow++;
                cout << s.b[pivotRow] << " ";
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
            for (int col = 1; col <= s.cols; col++)
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
                    cout << -s.A[pivotRow][freeColumn[freeIndex]] << " ";
                }
            }
            cout << ")T ";
        }
    }
    // cin.get();
    // cin.get();
    return FEASIBLE;
}
