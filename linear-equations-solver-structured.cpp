#include <cmath>
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
 * The work is split into four steps, each of which is one function:
 *     read the system        readSystem
 *     turn A into RREF       reduceToRref
 *     ask whether it is solvable
 *                            isInconsistent
 *     print the solution set printUniqueSolution / printGeneralSolution
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
 * The whole system, so that the functions below can work on it without
 * reaching for variables that live in main.
 */
struct System
{
    int rows = 0;                 // number of equations
    int cols = 0;                 // number of unknowns
    double A[MAXSIZE][MAXSIZE]{}; // coefficient matrix
    double b[MAXSIZE]{};          // constant vector
};

/*
 * The reduced row echelon form of a system, that is what the elimination made
 * out of it: the rank of the matrix, and the pivot column of every nonzero row.
 */
struct Rref
{
    int rank = 0;                 // number of nonzero rows
    int pivotOfRow[MAXSIZE]{};    // column of the leading entry of each of them
};

/*
 * The variables of the reduced system, split into the pivot ones, whose value
 * the constants give at once, and the free ones, which may be chosen at will.
 */
struct SolutionShape
{
    int freeCount = 0;            // number of free columns
    int freeColumn[MAXSIZE]{};    // the free columns themselves
    bool isFree[MAXSIZE]{};       // judge each column is or is not free
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

/*
 * Turn the system into reduced row echelon form and report its rank.
 *
 * This is the only step that changes the system, and it does so in two sweeps
 * of the same row operation:
 *     downward, to put a 1 in every pivot and a 0 below it, which is row
 *     echelon form;
 *     upward, to clear the entries above every pivot as well.
 *
 * The pivots are recorded while the first sweep finds them, so the second
 * sweep never has to search for the pivot column of a row again.
 */
Rref reduceToRref(System &s, double tolerance)
{
    Rref rref;

    // Forward elimination
    // Process one pivot after another. At the end of every round the pivot of
    // the current row is 1 and every entry below it is 0.
    for (int row = 1; row <= s.rows && rref.rank < s.cols; row++)
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
        rref.pivotOfRow[++rref.rank] = pivot; // record the pivot column

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

    // Back subtitution
    // Clear the entries above every pivot as well. Every pivot is already 1
    // and the pivot columns are known, so there is nothing to search for here.
    for (int row = rref.rank; row >= 1; row--)
    {
        // Eliminate other rows
        for (int target = 1; target <= row - 1; target++)
        {
            eliminate(s, row, target, rref.pivotOfRow[row]);
        }
    }

    return rref;
}

/*
 * Answer whether the system can be solved at all.
 *
 * In reduced row echelon form every row below the nonzero ones is a zero row,
 * so such a row means an equation "0 = b[row]"; the system admits a solution
 * iff none of those constants is nonzero.
 */
bool isInconsistent(const System &s, const Rref &rref, double tolerance)
{
    for (int row = rref.rank + 1; row <= s.rows; row++)
    {
        if (abs(s.b[row]) > tolerance)
        {
            return true;
        }
    }
    return false;
}

/*
 * Split the variables of the reduced system into pivot and free ones.
 *
 * A column is a pivot column iff some row has its first nonzero entry in it.
 * A free column carries no pivot, so its variable can be chosen freely and
 * becomes a parameter of the solution set.
 */
SolutionShape findFreeColumns(const System &s, const Rref &rref)
{
    SolutionShape shape;

    bool isPivot[MAXSIZE] = {}; // judge each column is or is not a pivot column
    for (int row = 1; row <= rref.rank; row++)
    {
        isPivot[rref.pivotOfRow[row]] = true;
    }

    for (int col = 1; col <= s.cols; col++)
    {
        if (!isPivot[col])
        {
            shape.isFree[col] = true;
            shape.freeCount++;
            shape.freeColumn[shape.freeCount] = col;
        }
        else
        {
            shape.isFree[col] = false;
        }
    }

    return shape;
}

/*
 * Print the only solution of a system whose every column carries a pivot, so
 * that the reduced system simply reads "x = the constant column".
 */
void printUniqueSolution(const System &s)
{
    cout << "The solution exists and is unique." << endl;
    cout << "Solution: ( ";
    for (int row = 1; row <= s.cols; row++)
    {
        cout << s.b[row] << " ";
    }
    cout << ")T";
}

/*
 * Print the whole solution set of a consistent system that has free variables.
 *
 * The particular solution takes every free variable to be 0, which leaves the
 * constant column for the pivot variables.
 *
 * Every basis vector takes one free variable to be 1 and the others to be 0;
 * the pivot entries of that vector are then minus the corresponding column of
 * the reduced matrix, because that is exactly what the equation of the pivot
 * row says about the free variable.
 */
void printGeneralSolution(const System &s, const SolutionShape &shape)
{
    cout << "There are infinitely many solutions." << endl;
    cout << "Solution: ( ";

    // Output constant part;
    // The pivot variables run through the rows in order, while the free
    // variables are set to zero
    int pivotRow = 0; // row that holds the pivot of the current column
    for (int col = 1; col <= s.cols; col++)
    {
        if (shape.isFree[col])
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
    for (int freeIndex = 1; freeIndex <= shape.freeCount; freeIndex++)
    {
        pivotRow = 0;
        cout << "+ x" << shape.freeColumn[freeIndex] << " ( ";
        for (int col = 1; col <= s.cols; col++)
        {
            if (shape.isFree[col])
            {
                if (col == shape.freeColumn[freeIndex])
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
                cout << -s.A[pivotRow][shape.freeColumn[freeIndex]] << " ";
            }
        }
        cout << ")T ";
    }
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

    // Reduce and judge
    double tolerance = zeroTolerance(s);
    Rref rref = reduceToRref(s, tolerance);
    if (isInconsistent(s, rref, tolerance))
    {
        cout << "There is no solution.";
        return INFEASIBLE;
    }

    // Output
    if (rref.rank == s.cols)
    {
        printUniqueSolution(s);
    }
    else
    {
        printGeneralSolution(s, findFreeColumns(s, rref));
    }
    // cin.get();
    // cin.get();
    return FEASIBLE;
}
