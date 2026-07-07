#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <regex>

int main() 
{

    const std::string inputFile = "horizons_results.txt";
    const std::string outputFile = "horizons_results_Saturn.txt";

    // Open input file for reading
    std::ifstream fin(inputFile);
    if (!fin.is_open()) 
    {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return 1;
    }

    // Open output file for writing
    std::ofstream fout(outputFile);
    if (!fout.is_open()) 
    {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        return 1;
    }


    fout << std::left << std::fixed << std::setw(20) << "\t\tX" << std::setw(20) << "\t\t Y \t\t" << "Z\n";

    // Use library Regex and function pattern to capture X, Y, and Z numeric values in a file
    // With this structue :
    //   X =-9.901283494896050E-01 Y = 1.278355489635163E-01 Z = 1.679920267418484E-04
    //
    // Explanation of the regex pattern:
    // Pattern used:
    //   X\s*=\s*([+-]?\d[.\deE+-]*)\s*Y\s*=\s*([+-]?\d[.\deE+-]*)\s*Z\s*=\s*([+-]?\d[.\deE+-]*)
    //
    // Each part captures a floating-point number (possibly in scientific notation) after the X =, Y =, and Z =.
    // Detail of each component:
    //
    //   X\s*=\s*       -> Matches the literal "X", followed by any number of spaces (\s*),
    //                    an equal sign "=", and again any number of optional spaces.
    //                  
    //
    //   ([+-]?\d[.\deE+-]*) -> This is the first capturing group that matches a float or scientific notation.
    //       [+-]?           -> It matches either a "+" or "-" sign.
    //       \d              -> Matches at least a single digit
    //       [.\deE+-]*      -> Matches zero or more of the following characters:
    //                      - digits (0–9),
    //                      - a decimal point ".",
    //                      - exponent marker "e" or "E",
    //                      - either '+' or '-' in the exponent.
    // 
    //  This allows matching both scientific notation like +2.34e+05 or -7.89E-03, and just decimal like 1.23 or -0.456
    //  For more information click on the link in annexe about regex library.                   
    //                      
    //                      
    //                      
    //
    //   Then we reuse the same pattern for Y and Z with their respective labels:
    //
    // 
    //   In total, this regex extracts three values from lines like:
    //   "X = -9.90E-01   Y = 1.27E-01   Z = 1.67E-04"
    //   Really useful to extract information from horizon file, without that it could have taken milion years to do that
    //   Resulting capturing groups:
    //   match[1] → value of X
    //   match[2] → value of Y
    //   match[3] → value of Z
    //
    std::regex pattern(R"(X\s*=\s*([+-]?\d[.\deE+-]*)\s*Y\s*=\s*([+-]?\d[.\deE+-]*)\s*Z\s*=\s*([+-]?\d[.\deE+-]*))");

    std::string line;
    while (std::getline(fin, line)) // read line till the end (empty line the ending)
    {
        // Trim leading spaces
        std::size_t startPos = line.find_first_not_of(" \t"); // The name speak for itself, it looks for the first char that is not a space or a tab
        if (startPos == std::string::npos) // ignore line if equal npos ( meaning that the line either is empty or whitespace-only)
        {
            continue;
        }

        // Clean space at the beginning
        std::string stripped = line.substr(startPos);

        // Check if this stripped line starts with "X =" 
        // if this the case return 0 and fullfil the condition
        if (stripped.rfind("X =", 0) == 0) 
        {
            // Attempt to match the regex to capture X, Y, Z
            std::smatch match;
            if (std::regex_search(stripped, match, pattern)) 
            {
                // match[1], match[2], match[3] are numeric strings
                if (match.size() == 4) 
                {
                    // Write them to the output file
                    fout << match[1] << " " 
                         << match[2] << " " 
                         << match[3] << "\n";
                }
            }
        }
    }

    fin.close();
    fout.close();

    std::cout << "Extraction complete. Results saved in " << outputFile << "\n";
    return 0;
}
