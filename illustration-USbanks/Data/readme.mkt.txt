Emir Malikov, Subal C. Kumbhakar, and Efthymios G. Tsionas, "A Cost 
System Approach to the Stochastic Directional Technology Distance 
Function with Undesirable Outputs: The Case of U.S. Banks in 2001-2010," 
Journal of Applied Econometrics, Vol. 31, No. 7, 2016, pp. 1407-1429.

The data file mkt2015-data.txt is a tab-delimited ASCII file in DOS
format. It is zipped in the file mkt-data.zip. Unix/Linux users should
use "unzip -a".

The file mkt2015-data.txt contains data for an unbalanced panel with
2,397 bank-year observations for 285 large U.S. commercial banks
operating in 2001-2010, whose total assets were in excess of one
billion dollars (in 2005 U.S. dollars) in the first three years of
observation. The data come from Call Reports available from the
Federal Reserve Bank of Chicago. For detailed description of the data
construction, please see Section 5 of the paper as well as the
Appendix of "Estimation of Banking Technology under Credit
Uncertainty" by Emir Malikov, Diego A. Restrepo-Tobon and Subal C.
Kumbhakar (2015, Empirical Economics, 49(1):185-211).

The list of included variables is as follows:
-----------------------------------------------------------------------
variable name   variable description
-----------------------------------------------------------------------
index           Bank ID
year            Year
y1              Consumer Loans, in real USD 1000
y2              Real Estate Loans, in real USD 1000
y3              Commercial & Industrial Loans, in real USD 1000
y4              Securities, in real USD 1000
y5              Off-Balance Sheet Activities Income, in real USD 1000
b               Total Reported Nonperforming Loans, in real USD 1000
x1              Labor, # of FT Employees
x2              Physical Capital (Fixed Assets), in real USD 1000
x3              Purchased Funds, in real USD 1000
x4              Interest-Bearing Transaction Accts, in real USD 1000
x5              Non-Transaction Accts, in real USD 1000
e               (Quasi-Fixed) Equity Capital, in real USD 1000
w1              Price of x1
w2              Price of x2
w3              Price of x3
w4              Price of x4
w5              Price of x5
TA              Total Assets
-----------------------------------------------------------------------

Data are sorted by index and year.

Nominal stock variables are deflated to 2005 USD using CPI Index (all
urban consumers).
