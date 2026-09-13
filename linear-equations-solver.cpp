#include <iostream>

using namespace std;

const double eps = 1e-7;

int main()
{
    // Define and input
    int m, n; // number of rows and columns of the coefficient matrix
    double A[102][102], b[102];
    cin >> m >> n;
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            cin >> A[i][j];
        }
        cin >> b[i];
    }

    // Forward elimination
    for (int i = 1; i <= m; i++)
    {
        int startj;   // the first nonzero column
        int maxr = i; // the row which has the biggest first element

        // Find nonzero column
        bool flag = false; // existence of nonzero column
        for (startj = i; startj <= n; startj++)
        {
            for (int ip = i; ip <= m; ip++) // ip means "iprime"
            {
                if (abs(A[ip][startj]) > abs(A[maxr][startj]))
                {
                    maxr = ip;
                }
            }
            if (abs(A[maxr][startj]) > eps)
            {
                flag = true;
                break;
            }
        }
        if (!flag)
        {
            break;
        }

        // Deal with the row with the biggest first element
        if (i != maxr)
        {
            swap(A[maxr], A[i]);
            swap(b[maxr], b[i]);
        }
        for (int j = startj + 1; j <= n; j++)
        {
            A[i][j] /= A[i][startj];
        }
        b[i] /= A[i][startj];
        A[i][startj] = 1;

        // Eliminate other rows
        for (int ip = i + 1; ip <= m; ip++)
        {
            for (int j = startj + 1; j <= n; j++)
            {
                A[ip][j] -= A[i][j] * A[ip][startj];
            }
            b[ip] -= b[i] * A[ip][startj];
            A[ip][startj] = 0;
        }
    }

    // Back subtitution
    for (int i = m; i >= 1; i--)
    {
        bool flag = false; // existence of nonzero element
        int nonzeroc;      // column of the nonzero element

        // Find nonzero column
        for (int j = 1; j <= n; j++)
        {
            if (abs(A[i][j]) > eps)
            {
                flag = true;
                nonzeroc = j;
                break;
            }
        }
        if (!flag)
        {
            if (abs(b[i]) > eps)
            {
                cout << "There is no solution.";
                return 0;
            }
            continue;
        }

        // Eliminate other rows
        for (int ip = 1; ip <= i - 1; ip++)
        {
            for (int j = nonzeroc + 1; j <= n; j++)
            {
                A[ip][j] -= A[i][j] * A[ip][nonzeroc];
            }
            b[ip] -= b[i] * A[ip][nonzeroc];
            A[ip][nonzeroc] = 0;
        }
    }

    // Judge the solution
    int nonzeror = 0; // number of nonzero rows
    for (int i = 1; i <= m; i++)
    {
        bool flag = false; // existence of nonzero element
        for (int j = 1; j <= n; j++)
        {
            if (abs(A[i][j]) > eps)
            {
                nonzeror++;
                flag = true;
                break;
            }
        }

        // Situation 1: no solution
        if (!flag)
        {
            if (abs(b[i]) > eps)
            {
                cout << "There is no solution.";
                return 0;
            }
        }
    }
    if (nonzeror == n)
    {
        // Situation 2: unique solution
        cout << "The solution exists and is unique." << endl;
        cout << "Solution: ( ";
        for (int i = 1; i <= n; i++)
        {
            cout << b[i] << " ";
        }
        cout << ")T";
    }
    else
    {
        // Situation 3: infinitely many solutions
        cout << "There are infinitely many solutions." << endl;

        // A column is a pivot column iff some row has its first nonzero entry in it
        bool ispivot[102];
        for (int j = 1; j <= n; j++)
        {
            ispivot[j] = false;
        }
        for (int i = 1; i <= nonzeror; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                if (abs(A[i][j]) > eps)
                {
                    ispivot[j] = true;
                    break;
                }
            }
        }

        // Find free columns (columns with free variables)
        int freec[102];   // free columns
        bool isfree[102]; // judge each column is or is not free
        int numfreec = 0; // number of free columns
        for (int j = 1; j <= n; j++)
        {
            if (!ispivot[j])
            {
                isfree[j] = true;
                numfreec++;
                freec[numfreec] = j;
            }
            else
            {
                isfree[j] = false;
            }
        }

        // Output
        cout << "Solution: ( ";

        // Output constant part;
        int count = 0;
        for (int c = 1; c <= n; c++) // c stands for "column"
        {
            if (isfree[c])
            {
                cout << "0 ";
            }
            else
            {
                count++;
                cout << b[count] << " ";
            }
        }
        cout << ")T ";

        // Output free part;
        for (int c = 1; c <= numfreec; c++)
        {
            count = 0;
            cout << "+ x" << freec[c] << " ( ";
            for (int j = 1; j <= n; j++)
            {
                if (isfree[j])
                {
                    if (j == freec[c])
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
                    count++;
                    cout << -A[count][freec[c]] << " ";
                }
            }
            cout << ")T ";
        }
    }
    // cin.get();
    // cin.get();
    return 0;
}