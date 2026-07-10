// fileLoader.js

// Универсальная функция для загрузки файлов
function loadFile(filePath, outputElementId) {
    console.log(`Starting fetch request for ${filePath}`);

    fetch(filePath)
        .then(response => {
            console.log(`Fetch response received for ${filePath}:`, response);
            if (!response.ok) {
                throw new Error(`Failed to load ${filePath}: ${response.status} ${response.statusText}`);
            }
            return response.text();
        })
        .then(text => {
            console.log(`File content loaded from ${filePath}:`, text);
            document.getElementById(outputElementId).innerHTML = marked.parse(text);
            document.getElementById(outputElementId).style.display = 'block'; // Показать содержимое файла
        })
        .catch(error => {
            document.getElementById(outputElementId).textContent = 'Error: ' + error.message;
            document.getElementById(outputElementId).style.display = 'block'; // Показать ошибку
            console.error(`Fetch error for ${filePath}:`, error); // Это выведет ошибку в консоль
        });
}
