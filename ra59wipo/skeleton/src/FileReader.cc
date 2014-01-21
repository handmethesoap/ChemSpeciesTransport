
#include <FileReader.hh>
#include <fstream>
#include <Debug.hh>
#include <iostream>

void FileReader::registerIntParameter(const std::string &key, int init)
{
    registerParameter<int>(key,init);
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
    registerParameter<real>(key,init);
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
    registerParameter<std::string>(key,init);
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
    registerParameter<std::string>(key,in);
}

void FileReader::setParameter(const std::string &key, real in)
{
    registerParameter<real>(key,in);
}

void FileReader::setParameter(const std::string &key, int in)
{
    registerParameter<int>(key,in);
}


bool FileReader::readFile(const std::string &name)
{
    std::ifstream paramFile(name, std::ifstream::in);
    if (!paramFile.is_open()) {
        return false;
    }
    std::string line;

    while(!paramFile.eof()) {
        getline (paramFile,line);
        if (!paramFile) continue;
        std::string delimiter ="#";
        std::string token = line.substr(0, line.find(delimiter));
        std::stringstream sst(token);
        std::string key, value;
        while (sst >> key) {
            if (sst >> value) {
                storage_[key]=value;
            } else {
                CHECK_MSG(false, "Parameter does not exist");
            }
        }
    }

    paramFile.close();
    return true;
}



void FileReader::printParameters() const
{
    for (auto iter=storage_.begin(); iter != storage_.end(); iter++) {
        std::cout << iter->first << "\t" << iter->second << std::endl;
    }
}
