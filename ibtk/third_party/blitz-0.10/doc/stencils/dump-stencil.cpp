#include <blitz/array.h>
#include <iostream>

BZ_USING_NAMESPACE(blitz)


BZ_DECLARE_STENCIL2(footprint,A,B)
    B = stenciltoapply(stencilargs);
BZ_END_STENCIL


template<class T_stencil>
void dumpStencil(const T_stencil& stencil)
{
    const int maxN = 5;

    Array<double,2> A(Range(-maxN,maxN),Range(-maxN,maxN)),
     B(Range(-maxN,maxN), Range(-maxN, maxN));

    A = 0.0;
    A(0,0) = 1;

    B = 0.0;

    applyStencil(stencil, A, B);

    B.reverseSelf(firstDim);
    B.reverseSelf(secondDim);

    using namespace blitz::tensor;

    int minrow = first(sum(B != 0.0, j));
    int maxrow = last(sum(B != 0.0, j));
    int mincol = first(sum(B(j,i) != 0.0, j));
    int maxcol = last(sum(B(j,i) != 0.0, j));

    // Code to dump HTML
		cout << "@html" << endl;
    cout << "<table cellpadding=2 rules=all>" << endl
         << "<tr align=right><td></td>";

    for (int j=mincol; j <= maxcol; ++j)
        cout << "<td>" << j << "</td>";

    cout << "</tr>";

    for (int i=minrow; i <= maxrow; ++i)
    {
        cout << "<tr align=right>";

        if (minrow != maxrow)
            cout << "<td>" << i << "</td>";
        else
            cout << "<td></td>";

        for (int j=mincol; j <= maxcol; ++j)
        {
            if (B(i,j) != 0.0)
            {
                if ((i == 0) && (j == 0))
                    cout << "<td bgcolor=\"#000060\">";
                else
                    cout << "<td bgcolor=\"#000000\">";
                cout << "<font color=\"#ffffff\">" << B(i,j) << "</font></td>";
            }
            else
                cout << "<td></td>";
        }

        cout << "</tr>" << endl;
    }
    cout << "</table>" << endl;
    cout << "@end html" << endl;

    // Code to dump TeX
		cout << "@tex" << endl;
		cout << "\\def\\tablevrule{\\vrule height 10pt depth 3pt}" << endl;
		cout << "\\vskip0.5\\baselineskip\\hskip\\itemindent" << endl;
		cout << "\\vbox{\\halign{\\hfill#\\quad&\\tablevrule#\\quad" << endl;
    for (int j=mincol; j <= maxcol; ++j)
        cout << "&\\hfill#\\quad";
    cout << "\\cr" << endl;
    cout << "&& ";
    for (int j=mincol; j < maxcol; ++j)
        cout << j << " & " ;
    cout << maxcol << " \\cr" << endl << "\\noalign{\\hrule}" << endl;
    for (int i=minrow; i <= maxrow; ++i)
    {
        if (minrow != maxrow)
            cout << i << " && ";
        else
            cout << " && ";
        for (int j=mincol; j <= maxcol; ++j)
        {
            if (B(i,j) != 0.0)
            {
                if ((i == 0) && (j == 0))
                    cout << "{\\bf " << B(i,j) << "} ";
                else
                    cout << " " << B(i,j) << " ";
            }
            if (j != maxcol)
                cout << " & ";
        }
        if (i != maxrow)
            cout << "\\cr" << endl;
    }
    cout <<  "\\cr}}" << endl;
    cout << "@end tex" << endl;

    // Code to dump !TeX && !HTML
		cout << "@ifnottex" << endl;
		cout << "@ifnothtml" << endl;
    cout << "@multitable {12345}";
    for (int j=mincol; j <= maxcol; ++j)
        cout << " {12345}";
    cout << endl << "@item @tab ";
    for (int j=mincol; j < maxcol; ++j)
        cout << j << " @tab " ;
    cout << maxcol << endl << "@item ";
    for (int i=minrow; i <= maxrow; ++i)
    {
        if (minrow != maxrow)
            cout << i << " @tab ";
        else
            cout << "@tab ";
        for (int j=mincol; j <= maxcol; ++j)
        {
            if (B(i,j) != 0.0)
            {
                if ((i == 0) && (j == 0))
                    cout << "@strong{" << B(i,j) << "} ";
                else
                    cout << " " << B(i,j) << " ";
            }
            if (j != maxcol)
                cout << "@tab ";
        }
        if (i != maxrow)
            cout << endl << "@item ";
    }
    cout << endl << "@end multitable" << endl;
		cout << "@end ifnothtml" << endl;
		cout << "@end ifnottex" << endl;

#if 0
    // Code to dump LaTeX
		cout << "latexcommand(" << endl;
    cout << "\\begin{tabular}{r|";
    for (int j=mincol; j <= maxcol; ++j)
        cout << "r";
    cout << "}" << endl;
    cout << " & ";
    for (int j=mincol; j < maxcol; ++j)
        cout << j << " & " ;
    cout << maxcol << " \\\\ \\hline" << endl;
    for (int i=minrow; i <= maxrow; ++i)
    {
        if (minrow != maxrow)
            cout << i << " & ";
        else
            cout << " & ";
        for (int j=mincol; j <= maxcol; ++j)
        {
            if (B(i,j) != 0.0)
            {
                if ((i == 0) && (j == 0))
                    cout << "{\\bf " << B(i,j) << "} ";
                else
                    cout << " " << B(i,j) << " ";
            }
            if (j != maxcol)
                cout << " & ";
        }
        if (i != maxrow)
            cout << " \\\\  " << endl;
    }
    cout << endl << "\\end{tabular}" << endl;
    cout << ")" << endl;
#endif 
}

int main()
{
    dumpStencil(footprint());
}

