import React from 'react';
import ashkanPic from './media/Ashkan.jpg';
import TarifPic from './media/Terrif.jpg';
import './styling/TwoCardsComponent.css';

const TwoCardsComponent = () => {
  const cards = [
    {
      id: 1,
      image: ashkanPic,
      title: 'Ashkan Nikfarjam',
      role: 'Lead Engineer',
      text: 'A Data Scinece studen, wiht passion for Data Analysis, Machine Learnin and AI.',
    },
    {
      id: 2,
      image: TarifPic,
      title: 'Tarif Khan',
      role: 'Front and Back end Engineer',
      text: 'I\'m studying computer science and passionate about programming, artificial intelligence and game design.',
    },
  ];

  return (
    <div className="cards-container">
      {cards.map((card) => (
        <div className="card" key={card.id}>
          <img className="card-image" src={card.image} alt={card.title} />
          <h2 className="card-title">{card.title}</h2>
          <h2 className="card-role">{card.role}</h2>
          <p className="card-text">{card.text}</p>
        </div>
      ))}
    </div>
  );
};

export default TwoCardsComponent;
