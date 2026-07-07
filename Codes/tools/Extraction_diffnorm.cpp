#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>

int main() 
{
    

    std::vector<std::string> planets = {"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"};
    

    std::ifstream infile("norm_diff.txt");
    if (!infile) 
    {
        std::cerr << "Error opening norm_diff.txt !" << std::endl;
        return 1;
    }
    
    // create a map that link planet to a output stream
    std::map<std::string, std::ofstream> outFiles;

    for (const auto& planet : planets) 
    {
        std::string filename = "norm_diff_Yoshida_" + planet + ".txt";
        outFiles[planet].open(filename);
        if (!outFiles[planet]) 
        {
            std::cerr << "Error opening file " << filename << std::endl;
            return 1;
        }
    }
    
    // Read the file norm_diff.txt line to line
    std::string line;

    while (std::getline(infile, line)) 
    {
        if (line.empty())
            {
                continue;
            }

        std::istringstream iss(line); // Split the line

    
        std::string col1, col2, col3, col4, col5; // Column
        // Extract each part and insert in variable ( char separator is a space)
        iss >> col1 >> col2 >> col3 >> col4 >> col5;

        

        if (!(iss)) // Not five column
        {
            
            continue;
        }
        
        // Search the current planet {planet}
        if (outFiles.find(col2) != outFiles.end()) 
        {
            outFiles[col2] << col1 << " " << col5 << "\n"; // Then write in the file mapped with the planet
        }
    }
    
    outFiles["Mercury"].close();
    outFiles["Venus"].close();
    outFiles["Earth"].close();
    outFiles["Mars"].close();
    outFiles["Jupiter"].close();
    
    infile.close();
    
    std::cout << "Extraction complete !" << std::endl;
    return 0;
}
