import React, { useEffect, useState } from 'react';

function App() {
  const [svgData, setSvgData] = useState(null);

  useEffect(() => {
      const fetchSvgData = async () => {
      try {
      const response = await fetch('http://localhost:8000');
      if (!response.ok) {
                  throw new Error('Network response was not ok');
                          }
      const text = await response.text();
      setSvgData(text);
      console.log( text )
      } catch (error) {
      console.error('Error fetching SVG data:', error);
      }
      };

      fetchSvgData();
      }, []);

  return (
      <div className="App">
      <h1>SVG Content:</h1>
      {svgData ? (
          <div dangerouslySetInnerHTML={{ __html: svgData }} />
          ) : (
            <p>Loading SVG data...</p>
            )}
      </div>
      );
}

export default App;
