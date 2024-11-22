import React from "react";
import { Bar } from "react-chartjs-2";

const BarChartConserved = ({ data }) => {
  const chartData = {
    labels: Array.from({ length: data.length }, (_, i) => `Position ${i + 1}`),
    datasets: [
      {
        label: "Conservation",
        data,
        backgroundColor: "rgba(75, 192, 192, 0.2)",
        borderColor: "rgba(75, 192, 192, 1)",
        borderWidth: 1,
      },
    ],
  };

  return (
    <div>
      <h2>Conserved Regions</h2>
      <Bar data={chartData} />
    </div>
  );
};

export default BarChartConserved;
