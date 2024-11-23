import React from "react";
import { Bar } from "react-chartjs-2";

const calculateBaseCounts = (sequences) => {
  return sequences.map((seq) => {
    const atCount = (seq.match(/A|T/gi) || []).length; // Count of A and T
    const gcCount = (seq.match(/G|C/gi) || []).length; // Count of G and C
    return { atCount, gcCount };
  });
};

const BarChartConserved = ({ sequences }) => {
  const counts = calculateBaseCounts(sequences);

  const chartData = {
    labels: sequences.map((_, index) => `Sequence ${index + 1}`),
    datasets: [
      {
        label: "AT Count",
        data: counts.map((item) => item.atCount),
        backgroundColor: "rgba(75,192,192,0.4)",
        borderColor: "rgba(75,192,192,1)",
        borderWidth: 1,
      },
      {
        label: "GC Count",
        data: counts.map((item) => item.gcCount),
        backgroundColor: "rgba(192,75,75,0.4)",
        borderColor: "rgba(192,75,75,1)",
        borderWidth: 1,
      },
    ],
  };

  return <Bar data={chartData} />;
};

export default BarChartConserved;
