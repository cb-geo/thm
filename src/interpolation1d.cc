#include <iomanip>
#include <iostream>
#include <vector>

double interpolate1d(std::vector<std::vector<double>>& Data, double x,
                     bool extrapolate) {

  std::vector<double> xData, yData;
  for (int i = 0; i < Data.size(); i++) {
    xData.push_back(Data[i][0]);
    yData.push_back(Data[i][1]);
  }

  int size = xData.size();

  for (int i = 0; i < size - 1; i++) {
    if (((x >= xData[i]) && (x <= xData[i + 1])) ||
        ((x <= xData[i]) && (x >= xData[i + 1]))) {
      double xL = xData[i], yL = yData[i], xR = xData[i + 1],
             yR = yData[i + 1];  // points on either side (unless beyond ends)
      double dydx = (yR - yL) / (xR - xL);  // gradient
      return yL + dydx * (x - xL);
      break;
    } else if ((x < xData[0]) && (x < xData[size - 1])) {
      (xData[0] < xData[size - 1] ? yData[0] : yData[size - 1]);
      break;
    } else if ((x > xData[0]) && (x > xData[size - 1])) {
      (xData[0] > xData[size - 1] ? yData[0] : yData[size - 1]);
      break;
    }
  }
}
