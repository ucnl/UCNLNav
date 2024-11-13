using System;

namespace UCNLNav.MML
{
    public class Matrix
    {
        #region Properties

        public static readonly double JG_ELIMINATE_EPS = 1E-5;

        int n, m;
        double[,] a;

        public int N { get { return n; } }
        public int M { get { return m; } }

        public bool IsSquare
        {
            get { return N == M; }
        }        

        #endregion

        #region Constructor

        public Matrix(int np, int mp)
        {
            if ((np <= 0) || (mp <= 0))
                throw new ArgumentOutOfRangeException("Dimensions should be greater than zero");

            n = np;
            m = mp;
            a = new double[n, m];
        }

        #endregion

        #region Methods

        public double this[int i, int j]
        {
            get { return a[i, j]; }
            set { a[i, j] = value; }
        }

        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            int m1n = m1.N;
            int m1m = m1.M;
            int m2m = m2.M;
            
            if (m2m != m1n)
                throw new ArgumentOutOfRangeException("m1.N should be equal to m2.M");

            Matrix result = new Matrix(m1n, m2m);

            for (int i = 0; i < m1n; i++)
            {
                for (int j = 0; j < m2m; j++)
                {
                    double sum = 0.0;

                    for (int k = 0; k < m1m; k++)
                    {
                        sum += m1.a[i, k] * m2.a[k, j];
                    }

                    result.a[i, j] = sum;
                }
            }

            return result;
        }

        public static Matrix Transpose(Matrix m)
        {
            Matrix result = new Matrix(m.M, m.N);

            for (int i = 0; i < m.N; i++)
            {
                for (int j = 0; j < m.M; j++)
                {
                    result.a[j, i] = m.a[i, j];
                }
            }

            return result;
        }

        /// <summary>
        /// Returns the augmented matrix as follows:
        ///     | a00 a01 a02 |                   | a00 a01 a02  1  0  0 |
        /// A = | a10 a11 a12 |, A.GetAugmented = | a10 a11 a12  0  1  0 |
        ///     | a20 a21 a22 |                   | a20 a21 a22  0  0  1 |
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static Matrix GetAugmented(Matrix m)
        {
            int rows_cnt = m.a.GetUpperBound(0) + 1;
            Matrix result = new Matrix(rows_cnt, rows_cnt * 2);

            for (int r = 0; r < rows_cnt; r++)
            {
                for (int c = 0; c < rows_cnt; c++)
                {
                    result.a[r, c] = m.a[r, c];
                }
                result.a[r, r + rows_cnt] = 1.0;
            }

            return result;
        }

        public static Matrix Inverse_JG(Matrix m)
        {
            if (!m.IsSquare)
                throw new ArgumentException("Matrix should be square");

            Matrix m_aug = Matrix.GetAugmented(m);

            int rows_cnt = m.a.GetUpperBound(0) + 1;
            int cols_cnt = rows_cnt * 2;

            for (int r = 0; r < rows_cnt; r++)
            {
                if (Math.Abs(m_aug.a[r, r]) < JG_ELIMINATE_EPS)
                {
                    for (int rr = r + 1; rr < rows_cnt; rr++)
                    {
                        if (Math.Abs(m_aug.a[rr, r]) > JG_ELIMINATE_EPS)
                        {
                            for (int c = 0; c < cols_cnt; c++)
                            {
                                double tmp = m_aug.a[r, c];
                                m_aug.a[r, c] = m_aug.a[rr, c];
                                m_aug.a[rr, c] = tmp;
                            }
                            break;
                        }
                    }
                }

                if (Math.Abs(m_aug.a[r, r]) > JG_ELIMINATE_EPS)
                {
                    for (int c = 0; c < cols_cnt; c++)
                    {
                        if (c != r)
                        {
                            m_aug.a[r, c] /= m_aug.a[r, r];
                        }
                    }
                    m_aug.a[r, r] = 1;

                    for (int rr = 0; rr < rows_cnt; rr++)
                    {
                        if (rr != r)
                        {
                            double factor = m_aug.a[rr, r] / m_aug.a[r, r];
                            for (int col = 0; col < cols_cnt; col++)
                            {
                                m_aug.a[rr, col] -= factor * m_aug.a[r, col];
                            }
                        }
                    }
                }
            }
            
            if (m_aug[rows_cnt - 1, rows_cnt - 1] == 0)
            {
                return null;
            }
            else
            {
                Matrix result = new Matrix(rows_cnt, rows_cnt);
                for (int r = 0; r < rows_cnt; r++)
                {
                    for (int c = 0; c < rows_cnt; c++)
                    {
                        result.a[r, c] = m_aug.a[r, c + rows_cnt];
                    }
                }

                return result;
            }
        }

        #endregion
    }
}
