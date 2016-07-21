#!/bin/sh
#
# This does some simple searches to look for common formatting issues,
# and a few potential bugs. As it uses regular expressions, it may have
# both false positives and false negatives.
#
# Usage: from top level of plasma directory:
# ./tools/checklist.sh
#
# See the ICL style guide
# https://bitbucket.org/icl/guides/wiki/icl_style_guide
#
# @author Mark Gates
# @date 2016-07-20

hg st -a -c -m | perl -pe 's/^[MAC] //' > files.txt

files=`cat files.txt | perl -pe 's/\n/ /'`
src=`grep -P '\.(c|h|cpp|hpp)\$' files.txt | perl -pe 's/\n/ /'`
hdr=`grep -P '\.(h|hpp)\$' files.txt | perl -pe 's/\n/ /'`

#echo "===== files\n$files"
#echo

#echo "===== src\n$src"
#echo

#echo "===== hdr\n$hdr"
#echo

function grep_src {
    echo "===== $1"
    #grep -c -P "$2" $src | grep -v ':0'
    grep -P "$2" $src
    echo
}

function grep_hdr {
    echo "===== $1"
    #grep -c -P "$2" $src | grep -v ':0'
    grep -P "$2" $hdr
    echo
}


# ---------------------------------------- format
echo "===== Windows returns. Use Unix newlines only."
grep -l -P '\r' $src
echo


# ---------------------------------------- indents
grep_src "Tabs. Use 4 spaces." '\t'

grep_src "Not 4 space indent. (This is really hard to search for.)" '^(    )* {1,3}(if|for|while|\{|\})'


# ---------------------------------------- spaces
grep_src "Trailing space." ' +$'

grep_src "Missing space before opening brace." '\S\{'

grep_src "Extra space before semicolon, comma." ' [,;]'

grep_src "Missing space after semicolon." ';[^ \n]'

# too many false positives like max(1,m) and A(i,j)
#grep_src "missing space after semicolon, comma" '[,;][^ \n]'

grep_src "Missing space after if, for, while." '^ +(if|for|while)\('

grep_src "Extra space inside if, for, while condition." '^ +(if|for|while) *\( '

# in (x ==y), can't check (x -y) because of unary - minus, nor (x *y) because of pointer declarations float *x.
echo "===== Relational (==, !=, <=, >=, <, >), assignment (+=, -=, *=, /=),"
echo "===== and boolean (&&, ||) operators should have spaces on both sides."
echo "===== Other operators may have spaces, which should be consistent, e.g., not (x +y)."
grep_src "Missing space around operators, e.g., (x==y)."  '([\w\[\]\(\)])(==|!=|<=|>=|<|>|\+=|-=|\*=|/=|&&|\|\|)([\w\[\]\(\)-])'
grep_src "Missing space after  operators, e.g., (x ==y)." '([\w\[\]\(\)]) +(==|!=|<=|>=|<|>|\+=|-=|\*=|/=|&&|\|\|\+|\/)([\w\[\]\(\)-])' | grep -v '#include <'
grep_src "Missing space before operators, e.g., (x== y)." '([\w\[\]\(\)])(==|!=|<=|>=|<|>|\+=|-=|\*=|/=|&&|\|\||\+|\-|\*|\/) +([\w\[\]\(\)-])'


# ---------------------------------------- newlines
grep_src "Cuddled curly braces: } else {. Add newline after {." '\} *else'

echo "===== Extra blank lines after {."
perl -n0777e 'if ( m/^(.*\{) *\n *\n/m ) { print "$ARGV: $1\n"; }' $src
echo

echo "===== Extra blank lines before }."
perl -n0777e 'if ( m/\n *\n( *\}.*)/ ) { print "$ARGV: $1\n"; }' $src
echo

echo "===== Missing newline at end of file (should be one)."
perl -n0777e 'if ( ! m/\n$/ ) { print "$ARGV\n"; }' $src
echo

echo "===== Extra blank lines at end of file (should be one)."
perl -n0777e 'print "$ARGV\n" if ( m/\n\n$/ )' $src
echo

echo "===== Missing blank line between functions (should be one)."
perl -n0777e 'if ( m/^\}\n\S/m ) { print "$ARGV\n"; }' $src
echo

echo "===== Extra blank lines between functions (should be one)."
perl -n0777e 'if ( m/^\}\n\n\n+\S/m ) { print "$ARGV\n"; }' $src
echo


# ---------------------------------------- potential bugs
echo "===== Headers not protected against multiple inclusion."
grep -P '^#ifndef \w+_([hH]_?|HPP)\b' -L $hdr
grep -P '^#define \w+_([hH]_?|HPP)\b' -L $hdr
echo

echo "===== Using REAL, COMPLEX, or PRECISION_[sdcz] without #define REAL, COMPLEX, or PRECISION."
grep -l -P 'defined\( *PRECISION_\w *\)'    $src | xargs grep -L -P '^ *#define PRECISION_\w'
grep -l -P 'ifdef +PRECISION_\w'            $src | xargs grep -L -P '^ *#define PRECISION_\w'
grep -l -P 'defined\( *(REAL|COMPLEX) *\)'  $src | xargs grep -L -P '^ *#define (REAL|COMPLEX)'
grep -l -P 'ifdef +(REAL|COMPLEX)'          $src | xargs grep -L -P '^ *#define (REAL|COMPLEX)'
echo


# ---------------------------------------- comments & docs
grep_src "C style /* ... */ comments (excluding Doxygen). Use C++ // style." '/\*[^*]'

grep_src "Using Fortran-style A**H or A**T. Please use A^H or A^T." '\*\*[TH].'
grep_src "Using Matlab-style A'. Please use A^H or A^T. (This is really hard to search for.)" " [A-Z]'"

echo "===== Using Matlab-style A'. (Second search.)"
grep -P "'" $src | grep -v -P "'[a-zA-Z=:+]'|'\\\0'"
echo

grep_src "Hermitian should be capitalized" "hermitian"