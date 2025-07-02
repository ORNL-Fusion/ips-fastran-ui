
import React, { useState, useEffect } from 'react';
import MonacoEditor from '@monaco-editor/react';

function MonacoFileEditor({ content }: { content: string }) {
  const [editorContent, setEditorContent] = useState(content);

  useEffect(() => {
    setEditorContent(content);
  }, [content]);

  const handleEditorChange = (newValue: string) => {
    setEditorContent(newValue);
  };

  return (
    <div style={{ height: '90vh' }}>
      <MonacoEditor
        width="100%"
        height="100%"
        language="python"
        theme="vs-dark"
        value={editorContent}
        options={{
          selectOnLineNumbers: true,
          renderWhitespace: "all",
        }}
      />
    </div>
  );
}

export default MonacoFileEditor;
