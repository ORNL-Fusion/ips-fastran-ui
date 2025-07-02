
import React, { useState } from 'react';

function FileExplorer({ onFileSelect }: { onFileSelect: (content: string, fileName: string) => void }) {
  const initialFiles = {
    'file1.py': 'print("Hello from file1.py")',
    'file2.py': 'print("Hello from file2.py")',
    'file3.py': 'print("Hello from file3.py")',
  };

  const [files, setFiles] = useState(Object.keys(initialFiles));
  const [fileContent, setFileContent] = useState<string>(initialFiles[files[0]]);
  const [fileName, setFileName] = useState<string>(files[0]);
  const [hoveredFile, setHoveredFile] = useState<string | null>(null); // Track hovered file

  const handleFileClick = (file: string) => {
    setFileContent(initialFiles[file]);
    setFileName(file);
    onFileSelect(initialFiles[file], file);
  };

  const handleDragStart = (e: React.DragEvent, index: number) => {
    e.dataTransfer.setData('text/plain', index.toString());
  };

  const handleDrop = (e: React.DragEvent, index: number) => {
    e.preventDefault();
    const draggedIndex = e.dataTransfer.getData('text/plain');
    const updatedFiles = [...files];
    const draggedFile = updatedFiles.splice(Number(draggedIndex), 1)[0];
    updatedFiles.splice(index, 0, draggedFile);
    setFiles(updatedFiles);
  };

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
  };

  return (
    <div className="file-explorer">
      <ul className="list-group">
        {files.map((file, index) => (
          <li
            key={index}
            className={`list-group-item bg-dark text-light file-item ${file === hoveredFile ? 'bg-secondary' : ''}`}
            draggable
            onDragStart={(e) => handleDragStart(e, index)}
            onDrop={(e) => handleDrop(e, index)}
            onDragOver={handleDragOver}
            onClick={() => handleFileClick(file)}
            onMouseEnter={() => setHoveredFile(file)} // Highlight on hover
            onMouseLeave={() => setHoveredFile(null)} // Remove highlight when hover ends
          >
            {file}
          </li>
        ))}
      </ul>
    </div>
  );
}

export default FileExplorer;
