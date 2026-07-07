#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

int main() 
{
    // Input and output file names
    const std::string file1 = "data_YoshidaMercury.txt";
    const std::string file2 = "horizons_results_Mercury.txt";
    const std::string diffFile = "differenceMercuryYoshida.txt";

    // Open the first input file
    std::ifstream fin1(file1);
    if (!fin1.is_open()) 
    {
        std::cerr << "Error opening reading file : " << file1 << std::endl;
        return 1;
    }

    // Open the second input file
    std::ifstream fin2(file2);
    if (!fin2.is_open()) 
    {
        std::cerr << "Error opening reading file : " << file2 << std::endl;
        return 1;
    }

    // Open the output file
    std::ofstream fout(diffFile);
    if (!fout.is_open()) 
    {
        std::cerr << "Error opening reading file : " << diffFile << std::endl;
        return 1;
    }

    // Skip the header line in both input files
    
    {
        std::string headerLine;
        if (!std::getline(fin1, headerLine)) 
        {
            std::cerr << "Error: " << file1 << " is empty or unreadable.\n";
            return 1;
        }
        if (!std::getline(fin2, headerLine)) 
        {
            std::cerr << "Error: " << file2 << " is empty or unreadable.\n";
            return 1;
        }
    }

    fout << std::left << std::fixed << "\tDX" << "\t\t\t\t" << "DY" << "\t\t\t\t" << "DZ\n";

    // Now read line by line from both files in parallel
    // Each line has three floating-point values: X, Y, Z
    std::string line1, line2;
    while (true) {
        if (!std::getline(fin1, line1)) 
        {
            break;
        }
        if (!std::getline(fin2, line2)) 
        {
            break;
        }

        // Parse the three values from line1
        double x1, y1, z1;
        {
            std::istringstream iss(line1);
            if (!(iss >> x1 >> y1 >> z1)) 
            {
                std::cerr << "Error parsing line in " << file1 << ": " << line1 << "\n";
                break;
            }
        }

        // Parse the three values from line2
        double x2, y2, z2;
        {
            std::istringstream iss(line2);
            if (!(iss >> x2 >> y2 >> z2)) 
            {
                std::cerr << "Error parsing line in " << file2 << ": " << line2 << "\n";
                break;
            }
        }

        // Compute differences
        double dx = x1 - x2;
        double dy = y1 - y2;
        double dz = z1 - z2;

        // Write the differences to the output file
        fout << dx << "\t\t" << dy << "\t\t" << dz << "\n";
    }

    fin1.close();
    fin2.close();
    fout.close();

    std::cout << "Differences saved to " << diffFile << "\n";
    return 0;
}
