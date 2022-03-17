const { sh, cli } = require("tasksfile");

function cleanAll() {
    sh("rm -rf build");
}

function createFieldSources() {
    sh("mkdir -p build");
    sh("npm install", {cwd: "depends/ffiasm"});
    sh("node ../depends/ffiasm/src/buildzqfield.js -q 21888242871839275222246405745257275088696311157297823662689037894645226208583 -n Fq", {cwd: "build"});
    sh("node ../depends/ffiasm/src/buildzqfield.js -q 21888242871839275222246405745257275088548364400416034343698204186575808495617 -n Fr", {cwd: "build"});
    
    if (process.platform === "darwin") {
        sh("nasm -fmacho64 --prefix _ fq.asm", {cwd: "build"});
    }  else if (process.platform === "linux") {
        sh("nasm -felf64 fq.asm", {cwd: "build"});
    } else throw("Unsupported platform");

    if (process.platform === "darwin") {
        sh("nasm -fmacho64 --prefix _ fr.asm", {cwd: "build"});
    }  else if (process.platform === "linux") {
        sh("nasm -felf64 fr.asm", {cwd: "build"});
    } else throw("Unsupported platform");
}

function buildEcdsa() {
    sh("g++" +
        " -I."+
        " -I../src"+
        " -I../depends/ffiasm/c"+
        " ../src/main.cpp"+
        " ../depends/ffiasm/c/misc.cpp"+
        " ../depends/ffiasm/c/naf.cpp"+
        " ../depends/ffiasm/c/splitparstr.cpp"+
        " ../depends/ffiasm/c/alt_bn128.cpp"+
        " fq.cpp"+
        " fq.o"+
        " fr.cpp"+
        " fr.o"+
        " -o ecdsa" +
        " -Wall -Wno-sign-compare"+
        " -O3 -fmax-errors=5 -std=c++17 -pthread -lgmp -lsodium", {cwd: "build", nopipe: true}
    );
}


cli({
    cleanAll,
    createFieldSources,
    buildEcdsa
});
