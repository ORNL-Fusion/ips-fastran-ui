
import React, { useState } from 'react';
import ActionList from './components/ActionList';
import FileExplorer from './components/FileExplorer';
import MonacoEditor from './components/MonacoEditor';

function App() {
  const [fileContent, setFileContent] = useState('');
  const [fileName, setFileName] = useState('');

  const handleFileSelect = (content: string, name: string) => {
    setFileContent(content);
    setFileName(name);
  };

  return (
    <div className="container-fluid d-flex vh-100">
      <div className="col-1 p-3 bg-dark">
        <ActionList />
      </div>
      <div className="col-2 p-3 bg-dark">
        <FileExplorer onFileSelect={handleFileSelect} />
      </div>
      <div className="col-9 p-3 bg-dark">
        <MonacoEditor content={fileContent} />
      </div>
    </div>
  );
}

export default App;
