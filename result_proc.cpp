#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


bool read_LOG_file(string file_name){
    string t; 
    string line;

    int s,w,r;
    int size_proc = 0;
    double area;
    double min_area = 0;
    int best_w;
    int best_r;
    bool new_porc = true;

    ifstream LR_file (file_name.c_str());
    if (LR_file.is_open())
    {
        while ( getline (LR_file,line) )
        {
            //cout << line << endl;
            istringstream l(line);
            l>>t;
            if (t == "config")
            {
                l>>t;
                l>>s;
                l>>t;
                l>>w;
                l>>t;
                l>>r;
                if (s != size_proc) //new RAM
                {
                    if (size_proc != 0)
                    {
                        cout << "For RAM size = " << size_proc << " best max width = " << best_w << " and ratio = " << best_r << " area = " << min_area << endl;
                    }
                    size_proc = s;
                    min_area = 0; //some max
                    new_porc = true;
                }
            }
            else if ( t== "Geometric")
            {
                //cout << line << endl;
                l>>t;
                l>>t;
                l>>t;
                l>>t;
                //l>>area;
                area = atof(t.c_str());
                //cout<< area << endl;
                //cout << min_area << endl;                
                if (area <= min_area || new_porc)
                {

                    best_w = w;
                    best_r = r;
                    min_area = area;
                    new_porc = false;
                }
            }
        }
	cout << "For RAM size = " << size_proc << " best max width = " << best_w << " and ratio = " << best_r << " area = " << min_area << endl;
        LR_file.close();
    }
    else 
    {
        cout << "unable to open file " << file_name << endl;
        return false;
    }
}



int main(int argc, char* argv[]) {
    string file_name;
    if (argc<2 )
    {
        cout << "check commandline!\n";
        return 0; 
    }
    file_name = argv[argc-1];

    read_LOG_file(file_name.c_str());

}
