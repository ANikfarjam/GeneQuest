import React from "react";
import { Bar } from "react-chartjs-2";
// import Chart from "chart.js/auto";

const BarChartConserved = ({ data }) => {
  const chartData = {
    labels: data.map((_, index) => `Position ${index}`),
    datasets: [
      {
        label: "Conservation Levels",
        data: data,
        backgroundColor: "rgba(75,192,192,0.4)",
        borderColor: "rgba(75,192,192,1)",
        borderWidth: 1,
      },
    ],
  };

  return <Bar data={chartData} />;
};

export default BarChartConserved;
