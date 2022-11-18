javac -d ./java/build ./java/*.java
jar cf ./build/Release/adaptagrams.jar -C ./java/build .
pause