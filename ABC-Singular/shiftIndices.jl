using Printf
powers = [
    0  0  0;
    0  0  1;
    1  0  0;
    0  1  0;
    0  0  2;
    1  0  1;
    0  1  1;
    1  1  0;
    2  0  0;
    0  2  0;
    0  0  3;
    1  0  2;
    0  1  2;
    1  1  1;
    2  0  1;
    0  2  1;
    2  1  0;
    1  2  0;
    3  0  0;
    0  3  0;
    0  0  4;
    1  0  3;
    0  1  3;
    1  1  2;
    2  0  2;
    0  2  2;
    2  1  1;
    1  2  1;
    2  2  0;
    3  0  1;
    0  3  1;
    3  1  0;
    1  3  0;
    4  0  0;
    0  4  0;
    0  0  5;
    1  0  4;
    0  1  4;
    1  1  3;
    2  0  3;
    0  2  3;
    2  1  2;
    1  2  2;
    2  2  1;
    3  0  2;
    0  3  2;
    3  1  1;
    1  3  1;
    3  2  0;
    2  3  0;
    4  0  1;
    0  4  1;
    4  1  0;
    1  4  0;
    5  0  0;
    0  5  0;
    0  0  6;
    1  0  5;
    0  1  5;
    1  1  4;
    2  0  4;
    0  2  4;
    2  1  3;
    1  2  3;
    2  2  2;
    3  0  3;
    0  3  3;
    3  1  2;
    1  3  2;
    3  2  1;
    2  3  1;
    3  3  0;
    4  0  2;
    0  4  2;
    4  1  1;
    1  4  1;
    4  2  0;
    2  4  0;
    5  0  1;
]
println(powers)
numberOfParameters = size(powers)[1]
println(numberOfParameters)

open("OTY7-shiftIndices.dat", "w") do outputFile::IOStream
    for i in 1:numberOfParameters
        @printf(outputFile, "%2d %2d %2d   %5d\n", powers[i, 1], powers[i, 2], powers[i, 3] + 7, i)
    end
end