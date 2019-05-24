
const std::map<std::string, int> create_measure_map()
{
  std::map<std::string, int> MEASURE;
  // {{"mic", 1}, {"mas", 2}, {"mev", 3}, {"mcn", 4}, {"tic", 5}, {"gmic", 6}}
  MEASURE["mic"]=1;
  MEASURE["mas"]=2;
  MEASURE["mev"]=3;
  MEASURE["mcn"]=4;
  MEASURE["tic"]=5;
  MEASURE["gmic"]=6;

  return MEASURE;
}

const std::map<std::string, int> create_est_map()
{
  std::map<std::string, int> EST;
  EST["mic_approx"]=0;
  EST["mic_e"]=1;
 
  return EST;
}