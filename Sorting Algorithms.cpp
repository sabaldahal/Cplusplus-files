#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <chrono>
#include <fstream>

typedef std::chrono::high_resolution_clock Clock;

#define ARRSIZE1 10
#define ARRSIZE2 100
#define ARRSIZE3 500
#define ARRSIZE4 5000
#define ARRSIZE5 25000
#define LINKEDLISTSIZE 50

using namespace std;

// bubble sort
void BubbleSort(int* array, int size) {
    for (int k = 0; k < (size-1); ++k) {
        int swapped = 0;
        for (int i = 0; i < (size - k - 1); ++i) {
            if (array[i] > array[i + 1]) {
                //swap items
                int temp = array[i];
                array[i] = array[i + 1];
                array[i + 1] = temp;
                swapped = 1;
            }
        }
        if (swapped == 0)
            break;
    }
}

//Insertion Sort
void InsertionSort(int* numbers, int size) {
    for (int i = 1; i < size; i++) {
        int j = i;
        while (j > 0 && numbers[j] < numbers[j - 1]) {
            // Swap numbers[j] and numbers [j - 1]
            int temp = numbers[j];
            numbers[j] = numbers[j - 1];
            numbers[j - 1] = temp;
            j--;
        }
    }
}

//Merge Sort
void Merge(int* numbers, int leftFirst, int leftLast, int rightLast) {
    int mergedSize = rightLast - leftFirst + 1;       // Size of merged partition
    int* mergedNumbers = new int[mergedSize]; // Dynamically allocates temporary
    // array for merged numbers
    int mergePos = 0;                         // Position to insert merged number
    int leftPos = leftFirst;                  // Initialize left partition position
    int rightPos = leftLast + 1;              // Initialize right partition position

    // Add smallest element from left or right partition to merged numbers
    while (leftPos <= leftLast && rightPos <= rightLast) {
        if (numbers[leftPos] <= numbers[rightPos]) {
            mergedNumbers[mergePos] = numbers[leftPos];
            leftPos++;
        }
        else {
            mergedNumbers[mergePos] = numbers[rightPos];
            rightPos++;
        }
        mergePos++;
    }

    // If left partition is not empty, add remaining elements to merged numbers
    while (leftPos <= leftLast) {
        mergedNumbers[mergePos] = numbers[leftPos];
        leftPos++;
        mergePos++;
    }

    // If right partition is not empty, add remaining elements to merged numbers
    while (rightPos <= rightLast) {
        mergedNumbers[mergePos] = numbers[rightPos];
        rightPos++;
        mergePos++;
    }

    // Copy merged numbers back to numbers
    for (mergePos = 0; mergePos < mergedSize; mergePos++) {
        numbers[leftFirst + mergePos] = mergedNumbers[mergePos];
    }

    // Free temporary array
    delete[] mergedNumbers;
}

void MergeSort(int* numbers, int startIndex, int endIndex) {
    if (startIndex < endIndex) {
        // Find the midpoint in the partition
        int mid = (startIndex + endIndex) / 2;

        // Recursively sort left and right partitions
        MergeSort(numbers, startIndex, mid);
        MergeSort(numbers, mid + 1, endIndex);

        // Merge left and right partition in sorted order
        Merge(numbers, startIndex, mid, endIndex);
    }
}

//Quick Sort
int Partition(int* numbers, int startIndex, int endIndex) {
    // Select the middle value as the pivot.
    int midpoint = startIndex + (endIndex - startIndex) / 2;
    int pivot = numbers[midpoint];

    // "low" and "high" start at the ends of the partition
    // and move toward each other.
    int low = startIndex;
    int high = endIndex;

    bool done = false;
    while (!done) {
        // Increment low while numbers[low] < pivot
        while (numbers[low] < pivot) {
            low = low + 1;
        }

        // Decrement high while pivot < numbers[high]
        while (pivot < numbers[high]) {
            high = high - 1;
        }

        // If low and high have crossed each other, the loop
        // is done. If not, the elements are swapped, low is
        // incremented and high is decremented.
        if (low >= high) {
            done = true;
        }
        else {
            int temp = numbers[low];
            numbers[low] = numbers[high];
            numbers[high] = temp;
            low = low + 1;
            high = high - 1;
        }
    }

    // "high" is the last index in the left partition.
    return high;
}

void Quicksort(int* numbers, int startIndex, int endIndex) {
    // Only sort if at least 2 elements exist
    if (endIndex <= startIndex) {
        return;
    }

    // Partition the array
    int high = Partition(numbers, startIndex, endIndex);

    // Recursively sort the left partition
    Quicksort(numbers, startIndex, high);

    // Recursively sort the right partition
    Quicksort(numbers, high + 1, endIndex);
}

//Heap Sort
void heapify(int* arr, int n, int i) {
    // Find largest among root, left child and right child
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && arr[left] > arr[largest])
        largest = left;

    if (right < n && arr[right] > arr[largest])
        largest = right;

    // Swap and continue heapifying if root is not largest
    if (largest != i) {
        swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

// main function to do heap sort
void heapSort(int* arr, int n) {
    // Build max heap
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // Heap sort
    for (int i = n - 1; i >= 0; i--) {
        swap(arr[0], arr[i]);

        // Heapify root element to get highest element at root again
        heapify(arr, i, 0);
    }
}


//Counting Sort
void countSort(int* array, int size) {

    int* output = new int[size+1];

    int max = array[0];

    // Find the largest element of the array
    for (int i = 1; i < size; i++) {
        if (array[i] > max)
            max = array[i];
    }
    int* count = new int[max + 1];
    // Initialize count array with all zeros.
    for (int i = 0; i <= max; ++i) {
        count[i] = 0;
    }

    // Store the count of each element
    for (int i = 0; i < size; i++) {
        count[array[i]]++;
    }

    // Store the cummulative count of each array
    for (int i = 1; i <= max; i++) {
        count[i] += count[i - 1];
    }

    // Find the index of each element of the original array in count array, and
    // place the elements in output array
    for (int i = size - 1; i >= 0; i--) {
        output[count[array[i]] - 1] = array[i];
        count[array[i]]--;
    }

    // Copy the sorted elements into original array
    for (int i = 0; i < size; i++) {
        array[i] = output[i];
    }
}

//Radix Sort
int RadixGetLength(int value) {
    if (value == 0) {
        return 1;
    }

    int digits = 0;
    while (value != 0) {
        digits++;
        value /= 10;
    }
    return digits;
}

// Returns the maximum length, in number of digits, out of all array elements
int RadixGetMaxLength(int* numbers, int numbersSize) {
    int maxDigits = 0;
    for (int i = 0; i < numbersSize; i++) {
        int digitCount = RadixGetLength(numbers[i]);
        if (digitCount > maxDigits) {
            maxDigits = digitCount;
        }
    }
    return maxDigits;
}

void RadixSort(int* numbers, int numbersSize) {
    vector<vector<int>*> buckets;
    for (int i = 0; i < 10; i++) {
        vector<int>* bucket = new vector<int>();
        buckets.push_back(bucket);
    }

    int copyBackIndex = 0;

    // Find the max length, in number of digits
    int maxDigits = RadixGetMaxLength(numbers, numbersSize);

    int pow10 = 1;
    for (int digitIndex = 0; digitIndex < maxDigits; digitIndex++) {
        // Put numbers into buckets
        for (int i = 0; i < numbersSize; i++) {
            int num = numbers[i];
            int bucketIndex = (abs(num) / pow10) % 10;
            buckets[bucketIndex]->push_back(num);
        }

        // Copy buckets back into numbers array
        copyBackIndex = 0;
        for (int i = 0; i < 10; i++) {
            vector<int>& bucket = *buckets[i];
            for (int j = 0; j < bucket.size(); j++) {
                numbers[copyBackIndex] = bucket[j];
                copyBackIndex++;
            }
            bucket.clear();
        }

        pow10 *= 10;
    }

    vector<int> negatives;
    vector<int> nonNegatives;
    for (int i = 0; i < numbersSize; i++) {
        int num = numbers[i];
        if (num < 0) {
            negatives.push_back(num);
        }
        else {
            nonNegatives.push_back(num);
        }
    }

    // Copy sorted content to array - negatives in reverse, then non-negatives
    copyBackIndex = 0;
    for (int i = negatives.size() - 1; i >= 0; i--) {
        numbers[copyBackIndex] = negatives[i];
        copyBackIndex++;
    }
    for (int i = 0; i < nonNegatives.size(); i++) {
        numbers[copyBackIndex] = nonNegatives[i];
        copyBackIndex++;
    }

    // Free each dynamically allocated bucket
    for (int i = 0; i < 10; i++) {
        delete buckets[i];
    }
}


//print array
void printArray(int* arr, int n) {
    for (int i = 0; i < n; ++i)
        cout << arr[i] << " ";
    cout << "\n";
}

/*
* Student class that stores student's information
*/
class Student {  
public:
    string firstName;
    string lastName;
    int mNumber;
    Student() {};
    Student(string firstName, string lastName, int mNumber) {
        this->firstName = firstName;
        this->lastName = lastName;
        this->mNumber = mNumber;
    }
};
/*
* Node class to point to students data and next and prev data in the linked list
*/
class Node {
public:
    Student* student;
    Node* prev;
    Node* next;


    Node() {};
    Node(Student* student) {
        this->student = student;
    }

};

/*
* Linked List class to manage every nodes pointing to students data
*/
class linkedList {
private:
    Node* head;
    Node* tail;
    int count{ 0 };
    char sortUsing;

public:
    linkedList(){}
    //adds item to the list
    int getCount() {
        return this->count;
    }
    void generateList() {
        int nlength; //length of either first or last name
        int mNumber; //random mnumber of 8 digits
        string f_name = "";
        string l_name = "";
        int m_num;
        char randomChar;
        for (int i = 0; i < LINKEDLISTSIZE; i++) {
            f_name = "";
            l_name = "";
            nlength = 3 + (rand() % 6); //random int for length of first name;
            //generating first name
            for (int j = 0; j < nlength; j++) {
                randomChar = 'a' + (rand() % 26);
                f_name += randomChar;
            }
            nlength = 3 + (rand() % 6); //random int for length of last name;
            //generating last name
            for (int k = 0; k < nlength; k++) {
                randomChar = 'a' + (rand() % 26);
                l_name += randomChar;
            }
            mNumber = 10000000 + (rand() % 89999999); //random number of 8 digits for mNumber

            //saving the student's info to the list
            Student* newStudent = new Student(f_name, l_name, mNumber);
            Node* newNode = new Node(newStudent);
            this->addItem(newNode);
        }

    }
    void addItem(Node* stu) {
        if (!this->head) {
            this->head = stu;
            this->tail = stu;
            this->count = 1;
            return;
        }
        this->tail->next = stu;
        stu->prev = this->tail;
        this->tail = stu;
        this->count++;
    }
    //compare two nodes, return true if the first node is greater
    bool isGreater(Node* item, Node* rhs) {

        int index = 0;
        bool done = false;
        int minSize = item->student->firstName.length();
        //sort using first Name
        if (this->sortUsing == 'f') {
            //calculate the smallest sized string and record the length
            if (minSize > rhs->student->firstName.length()) minSize = rhs->student->firstName.length();
            while (!done) {
                //handling index out of range
                if (index == minSize) return false;
                //compare and return true if this->student is smaller
                if (item->student->firstName.at(index) > rhs->student->firstName.at(index)) {
                    done = true;
                    return true;
                }
                //if the characters are same loop again
                else if (item->student->firstName.at(index) == rhs->student->firstName.at(index)) {
                    index++;
                }
                //return false if rhs is smaller
                else {
                    done = true;
                    return false;
                }
            }
        }
        //sort using last Name
        else if (this->sortUsing == 'l') {
            //calculate the smallest sized string and record the length
            minSize = item->student->lastName.length();
            if (minSize > rhs->student->lastName.length()) minSize = rhs->student->lastName.length();
            while (!done) {
                //handling index out of range
                if (index == (minSize)) return false;
                //compare and return true if this->student is smaller
                if (item->student->lastName.at(index) > rhs->student->lastName.at(index)) {
                    done = true;
                    return true;
                }
                //if the characters are same loop again
                else if (item->student->lastName.at(index) == rhs->student->lastName.at(index)) {
                    index++;
                }
                //return false if rhs is smaller
                else {
                    done = true;
                    return false;
                }
            }
        }
        //sort using MNumber
        else if (this->sortUsing == 'm') {
            if (item->student->mNumber > rhs->student->mNumber) return true;
            return false;
        }
    }
    //compare two nodes, return true if the first node is smaller
    bool isSmaller(Node* item, Node* rhs) {
        int index = 0;
        bool done = false;
        int minSize = item->student->firstName.length();
        //sort using first Name
        if (sortUsing == 'f') {
            //calculate the smallest sized string and record the length
            if (minSize > rhs->student->firstName.length()) minSize = rhs->student->firstName.length();
            while (!done) {
                //handling index out of range
                if (index == (minSize)) return false;
                //compare and return true if this->student is smaller
                if (item->student->firstName.at(index) < rhs->student->firstName.at(index)) {
                    done = true;
                    return true;
                }
                //if the characters are same loop again
                else if (item->student->firstName.at(index) == rhs->student->firstName.at(index)) {
                    index++;
                }
                //return false if rhs is smaller
                else {
                    done = true;
                    return false;
                }
            }
        }
        //sort using last Name
        else if (sortUsing == 'l') {
            //calculate the smallest sized string and record the length
            minSize = item->student->lastName.length();
            if (minSize > rhs->student->lastName.length()) minSize = rhs->student->lastName.length();
            while (!done) {
                //handling index out of range
                if (index == (minSize)) return false;
                //compare and return true if this->student is smaller
                if (item->student->lastName.at(index) < rhs->student->lastName.at(index)) {
                    done = true;
                    return true;
                }
                //if the characters are same loop again
                else if (item->student->lastName.at(index) == rhs->student->lastName.at(index)) {
                    index++;
                }
                //return false if rhs is smaller
                else {
                    done = true;
                    return false;
                }
            }
        }
        //sort using MNumber
        else if (sortUsing == 'm') {
            if (item->student->mNumber < rhs->student->mNumber) return true;
            return false;
        }
    }
    //print the entire list
    void printList() {
        Node* curr = this->head;
        int counter = 1;
        while (curr) {
            cout << counter << " ";
            cout << "Info: " << curr->student->firstName << " ";
            cout << curr-> student->lastName << " ";
            cout << curr->student->mNumber << endl;
            counter++;
            curr = curr->next;
        }
    }
    //swap two items in the list
    //swap consecutive items
    void clearList() {
        Node* curr = this->head;
        while (curr) {
            this->head = curr->next;
            delete curr;
            curr = this->head;
        }
        this->head = nullptr;
        this->tail = nullptr;
        this->count = 0;
    }
    void swapItems(Node* item1, Node* item2) {  
        if (item1 == this->head) this->head = item2;
        if (item2 == this->tail) this->tail = item1;
        item1->next = item2->next;
        if (item2->next) item2->next->prev = item1;
        item2->next = item1;
        item2->prev = item1->prev;
        if (item1->prev) item1->prev->next = item2;
        item1->prev = item2;     
    }
    //swap items that are not consecutive
    void swapItems2(Node* item1, Node* item2) {
        Node* temp1 = item1->next;
        Node* temp2 = item2->prev;
        if (item1 == this->head) this->head = item2;
        if (item2 == this->tail) this->tail = item1;
        item1->next = item2->next;
        if (item2->next) item2->next->prev = item1;
        item2->next = temp1;
        item2->prev = item1->prev;
        if (item1->prev) item1->prev->next = item2;
        item1->prev = temp2;
        temp1->prev = item2;
        temp2->next = item1;
    }
    //sort the list using Bubble Sort algorithm
    void BubbleSort(char sortBasis = 'f', char sortOrder = 'a') {
        if (!(sortBasis == 'f' || sortBasis == 'l' || sortBasis == 'm')) {
            cout << "Invalid Input";
            return;
        }
        if (!(sortOrder == 'a' || sortOrder == 'd')) {
            cout << "Invalid Sorting Order";
            return;
        }
        this->sortUsing = sortBasis;
        Node* curr = this->head;
        //ascending order
        if (sortOrder == 'a') {
            for (int k = 0; k < (this->count - 1); ++k) {
                curr = this->head;
                int swapped = 0;
                for (int i = 0; i < (this->count - k - 1); ++i) {
                    //recalculating the current node
                    curr = this->head;
                    for (int t = 0; t < i; t++) {
                        curr = curr->next;
                    }
                    //ascending order
                    if (isGreater(curr, curr->next)) {
                        //swap items
                        this->swapItems(curr, curr->next);
                        swapped = 1;
                    }
                }
                if (swapped == 0)
                    break;
            }
        }
        //descending order
        else {
            for (int k = 0; k < (this->count - 1); ++k) {
                curr = this->head;
                int swapped = 0;
                for (int i = 0; i < (this->count - k - 1); ++i) {
                    //recalculating the current node
                    curr = this->head;
                    for (int t = 0; t < i; t++) {
                        curr = curr->next;
                    }
                    //ascending order
                    if (isSmaller(curr, curr->next)) {
                        //swap items
                        this->swapItems(curr, curr->next);
                        swapped = 1;
                    }

                }
                if (swapped == 0)
                    break;
            }
        }



    }
    //sort the list using Insertion Sort algorithm
    void InsertionSort(char sortBasis = 'f', char sortOrder = 'a') {
        //asserting correct values
        if (!(sortBasis == 'f' || sortBasis == 'l' || sortBasis == 'm')) {
            cout << "Invalid Input";
            return;
        }
        if (!(sortOrder == 'a' || sortOrder == 'd')) {
            cout << "Invalid Sorting Order";
            return;
        }
        this->sortUsing = sortBasis;
        Node* curr = this->head->next;
        //ascending order
        if (sortOrder == 'a') {
            for (int i = 1; i < this->count; i++) {
                int j = i;
                //recalculating the current node
                curr = this->head;
                for (int t = 0; t < j; t++) {
                    curr = curr->next;
                }
                while (j > 0 && isSmaller(curr, curr->prev)) {
                    // swap items
                    this->swapItems(curr->prev, curr);
                    j--;
                    //recalculating the current node
                    curr = this->head;
                    for (int t = 0; t < j; t++) {
                        curr = curr->next;
                    }
                }
            }
        }
        //descending order
        else {
            for (int i = 1; i < this->count; i++) {
                int j = i;
                //recalculating the current node
                curr = this->head;
                for (int t = 0; t < j; t++) {
                    curr = curr->next;
                }
                while (j > 0 && isGreater(curr, curr->prev)) {
                    // swap items
                    this->swapItems(curr->prev, curr);
                    j--;
                    //recalculating the current node
                    curr = this->head;
                    for (int t = 0; t < j; t++) {
                        curr = curr->next;
                    }
                }
            }
        }

    }
    //sort the list using Quick Sort algorithm
//Quick Sort
    int Partition(int startIndex, int endIndex, char sortBasis, char sortOrder) {
        // Select the middle value as the pivot.
        int midpoint = startIndex + (endIndex - startIndex) / 2;
        Node* pivot = this->head;
        for (int i = 0; i < midpoint; i++) {
            pivot = pivot->next;
        }

        // "low" and "high" start at the ends of the partition
        // and move toward each other.
        int low = startIndex;
        int high = endIndex;


        bool done = false;
        while (!done) {
            //updating low Node with the respective address
            Node* lowNode = this->head;
            for (int itr = 0; itr < low; itr++) {
                lowNode = lowNode->next;
            }
            //updating high node with the respective address
            Node* highNode = this->head;
            for (int itr = 0; itr < high; itr++) {
                highNode = highNode->next;
            }
            //ascending order
            if (sortOrder == 'a') {
                // Increment low while lownode < pivot
                while (isSmaller(lowNode, pivot)) {
                    low = low + 1;
                    //updating low Node with the respective address
                    lowNode = this->head;
                    for (int itr = 0; itr < low; itr++) {
                        lowNode = lowNode->next;
                    }
                }

                // Decrement high while pivot < highnode
                while (isSmaller(pivot, highNode)) {
                    high = high - 1;
                    //updating high node with the respective address
                    highNode = this->head;
                    for (int itr = 0; itr < high; itr++) {
                        highNode = highNode->next;
                    }
                }
            }
            //descending order
            else {
                // Increment low while lownode > pivot
                while (isGreater(lowNode, pivot)) {
                    low = low + 1;
                    //updating low Node with the respective address
                    lowNode = this->head;
                    for (int itr = 0; itr < low; itr++) {
                        lowNode = lowNode->next;
                    }
                }

                // Decrement high while pivot > highnode
                while (isGreater(pivot, highNode)) {
                    high = high - 1;
                    //updating high node with the respective address
                    highNode = this->head;
                    for (int itr = 0; itr < high; itr++) {
                        highNode = highNode->next;
                    }
                }
            }



            // If low and high have crossed each other, the loop
            // is done. If not, the elements are swapped, low is
            // incremented and high is decremented.
            if (low >= high) {
                done = true;
            }
            else {
                //swap items
                //updating low Node with the respective address
                lowNode = this->head;
                for (int itr = 0; itr < low; itr++) {
                    lowNode = lowNode->next;
                }
                //updating high node with the respective address
                highNode = this->head;
                for (int itr = 0; itr < high; itr++) {
                    highNode = highNode->next;
                }
                int difference = abs(high - low);
                if (difference == 1) {
                    this->swapItems(lowNode, highNode);
                }
                else {
                    this->swapItems2(lowNode, highNode);
                }
                low = low + 1;
                high = high - 1;
            }
        }

        // "high" is the last index in the left partition.
        return high;
    }

    void Quicksort(int startIndex, int endIndex, char sortBasis = 'f', char sortOrder = 'a') {
        //asserting correct values
        if (!(sortBasis == 'f' || sortBasis == 'l' || sortBasis == 'm')) {
            cout << "Invalid Input";
            return;
        }
        if (!(sortOrder == 'a' || sortOrder == 'd')) {
            cout << "Invalid Sorting Order";
            return;
        }
        this->sortUsing = sortBasis;
        // Only sort if at least 2 elements exist
        if (endIndex <= startIndex) {
            return;
        }

        // Partition the array
        int high = Partition(startIndex, endIndex, sortBasis, sortOrder);

        // Recursively sort the left partition
        Quicksort(startIndex, high, sortBasis, sortOrder);

        // Recursively sort the right partition
        Quicksort(high + 1, endIndex, sortBasis, sortOrder);
    }
    
};

void testLinkedList() {
    //linked list test;
    linkedList newList;
    newList.generateList();
    cout << "Here is a list of the students" << endl;
    newList.printList();
    
    cout << endl;
    cout << "Choose one of the sorting algorithms: " << endl;
    cout << "1. Bubble Sort " << endl;
    cout << "2. Insertion Sort " << endl;
    cout << "3. Quick Sort " << endl;
    cout << endl;
    cout << "Please enter the number to choose one of the option: (1/2/3)? ";
    int option;
    char sortbasis;
    char sortorder;
    cin >> option;
    cout << endl;
    cout << "Please specify how you would like to sort the list: " << endl;
    cout << "Type (f) to sort by First Name " << endl;
    cout << "Type (l) to sort by Last Name " << endl;
    cout << "Type (m) to sort by MNumber " << endl;
    cout << endl;
    cout << "Enter your answer here: ";
    cin >> sortbasis;
    cout << endl;
    cout << "Please specify how you would like to sort the list: " << endl;
    cout << "Type (a) to sort in ascending order " << endl;
    cout << "Type (d) to sort in descending order " << endl;
    cout << endl;
    cout << "Enter your answer here: ";
    cin >> sortorder;
    cout << endl;
    
    //select one of the algorithms
    if (option == 1) {
        newList.BubbleSort(sortbasis, sortorder);
    }
    else if (option == 2) {
        newList.InsertionSort(sortbasis, sortorder);
    }
    else if (option == 3) {
        newList.Quicksort(0, newList.getCount() - 1, sortbasis, sortorder);
    }
    else {
        return;
    }
    newList.printList();
}

//test for all Sorting algorithms
int testBubbleSort(int* arr, int size) {
    //making a copy of the array
    int* unsortedArray = new int[size];
    for (int i = 0; i < size; i++) {
        unsortedArray[i] = arr[i];
    }
    //start calculating time
    auto t1 = Clock::now();
    BubbleSort(unsortedArray, size);
    auto t2 = Clock::now();
    delete[] unsortedArray;
    return chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
}
int testInsertionSort(int* arr, int size) {
    //making a copy of the array
    int* unsortedArray = new int[size];
    for (int i = 0; i < size; i++) {
        unsortedArray[i] = arr[i];
    }
    //start calculating time
    auto t1 = Clock::now();
    InsertionSort(unsortedArray, size);
    auto t2 = Clock::now();
    delete[] unsortedArray;
    return chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
}
int testMergeSort(int* arr, int size) {
    //making a copy of the array
    int* unsortedArray = new int[size];
    for (int i = 0; i < size; i++) {
        unsortedArray[i] = arr[i];
    }
    //start calculating time
    auto t1 = Clock::now();
    MergeSort(unsortedArray, 0, size-1);
    auto t2 = Clock::now();
    delete[] unsortedArray;
    return chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
}
int testQuickSort(int* arr, int size) {
    //making a copy of the array
    int* unsortedArray = new int[size];
    for (int i = 0; i < size; i++) {
        unsortedArray[i] = arr[i];
    }
    //start calculating time
    auto t1 = Clock::now();
    Quicksort(unsortedArray, 0 , size-1);
    auto t2 = Clock::now();
    delete[] unsortedArray;
    return chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
}
int testHeapSort(int* arr, int size) {
    //making a copy of the array
    int* unsortedArray = new int[size];
    for (int i = 0; i < size; i++) {
        unsortedArray[i] = arr[i];
    }
    //start calculating time
    auto t1 = Clock::now();
    heapSort(unsortedArray, size);
    auto t2 = Clock::now();
    delete[] unsortedArray;
    return chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
}
int testCountSort(int* arr, int size) {
    //making a copy of the array
    int* unsortedArray = new int[size];
    for (int i = 0; i < size; i++) {
        unsortedArray[i] = arr[i];
    }
    //start calculating time
    auto t1 = Clock::now();
    countSort(unsortedArray, size);
    auto t2 = Clock::now();
    delete[] unsortedArray;
    return chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
}
int testRadixSort(int* arr, int size) {
    //making a copy of the array
    int* unsortedArray = new int[size];
    for (int i = 0; i < size; i++) {
        unsortedArray[i] = arr[i];
    }
    //start calculating time
    auto t1 = Clock::now();
    RadixSort(unsortedArray, size);
    auto t2 = Clock::now();
    delete[] unsortedArray;
    return chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
}
int main()
{
    //run linked list implementation
    testLinkedList();
    //check if the user wants to test the run time complexity of the sorting algorithms
    char answer;
    cout << "Would you like to test the time for sorting algorithms? (y/n) ";
    cin >> answer;
    if (answer != 'y') return 0;

    //test time taken for every sorting algorithm
    int timetaken;
    ofstream outputFile;
    outputFile.open("./AllData.csv", ios_base::app);
   
	//create arrays and populate it
    //these arrays will remain unchanged;
    int* arr1 = new int[ARRSIZE1];
    int* arr2 = new int[ARRSIZE2];
    int* arr3 = new int[ARRSIZE3];
    int* arr4 = new int[ARRSIZE4];
    int* arr5 = new int[ARRSIZE5];
    
    for (int i = 0; i < ARRSIZE1; i++) {
        arr1[i] = rand() % (2 * ARRSIZE1);
    }
    for (int i = 0; i < ARRSIZE2; i++) {
        arr2[i] = rand() % (2 * ARRSIZE2);
    }
    for (int i = 0; i < ARRSIZE3; i++) {
        arr3[i] = rand() % (2 * ARRSIZE3);
    }
    for (int i = 0; i < ARRSIZE4; i++) {
        arr4[i] = rand() % (2 * ARRSIZE4);
    }
    for (int i = 0; i < ARRSIZE5; i++) {
        arr5[i] = rand() % (2 * ARRSIZE5);
    }
    //test each of the algorithm using all array size
    outputFile << "\n";
    //testing Bubble Sort and measuring time in ns
    outputFile << "BubbleSort,";
    timetaken = testBubbleSort(arr1, ARRSIZE1);
    outputFile << timetaken << ",";
    timetaken = testBubbleSort(arr2, ARRSIZE2);
    outputFile << timetaken << ",";
    timetaken = testBubbleSort(arr3, ARRSIZE3);
    outputFile << timetaken << ",";
    timetaken = testBubbleSort(arr4, ARRSIZE4);
    outputFile << timetaken << ",";
    timetaken = testBubbleSort(arr5, ARRSIZE5);
    outputFile << timetaken << "\n";

    //testing Insertion Sort
    outputFile << "InsertionSort,";
    timetaken = testInsertionSort(arr1, ARRSIZE1);
    outputFile << timetaken << ",";
    timetaken = testInsertionSort(arr2, ARRSIZE2);
    outputFile << timetaken << ",";
    timetaken = testInsertionSort(arr3, ARRSIZE3);
    outputFile << timetaken << ",";
    timetaken = testInsertionSort(arr4, ARRSIZE4);
    outputFile << timetaken << ",";
    timetaken = testInsertionSort(arr5, ARRSIZE5);
    outputFile << timetaken << "\n";

    //testing Merge Sort
    outputFile << "MergeSort,";
    timetaken = testMergeSort(arr1, ARRSIZE1);
    outputFile << timetaken << ",";
    timetaken = testMergeSort(arr2, ARRSIZE2);
    outputFile << timetaken << ",";
    timetaken = testMergeSort(arr3, ARRSIZE3);
    outputFile << timetaken << ",";
    timetaken = testMergeSort(arr4, ARRSIZE4);
    outputFile << timetaken << ",";
    timetaken = testMergeSort(arr5, ARRSIZE5);
    outputFile << timetaken << "\n";

    //testing Quick Sort
    outputFile << "QuickSort,";
    timetaken = testQuickSort(arr1, ARRSIZE1);
    outputFile << timetaken << ",";
    timetaken = testQuickSort(arr2, ARRSIZE2);
    outputFile << timetaken << ",";
    timetaken = testQuickSort(arr3, ARRSIZE3);
    outputFile << timetaken << ",";
    timetaken = testQuickSort(arr4, ARRSIZE4);
    outputFile << timetaken << ",";
    timetaken = testQuickSort(arr5, ARRSIZE5);
    outputFile << timetaken << "\n";

    //testing Heap Sort
    outputFile << "HeapSort,";
    timetaken = testHeapSort(arr1, ARRSIZE1);
    outputFile << timetaken << ",";
    timetaken = testHeapSort(arr2, ARRSIZE2);
    outputFile << timetaken << ",";
    timetaken = testHeapSort(arr3, ARRSIZE3);
    outputFile << timetaken << ",";
    timetaken = testHeapSort(arr4, ARRSIZE4);
    outputFile << timetaken << ",";
    timetaken = testHeapSort(arr5, ARRSIZE5);
    outputFile << timetaken << "\n";

    //testing Count Sort
    outputFile << "CountSort,";
    timetaken = testCountSort(arr1, ARRSIZE1);
    outputFile << timetaken << ",";
    timetaken = testCountSort(arr2, ARRSIZE2);
    outputFile << timetaken << ",";
    timetaken = testCountSort(arr3, ARRSIZE3);
    outputFile << timetaken << ",";
    timetaken = testCountSort(arr4, ARRSIZE4);
    outputFile << timetaken << ",";
    timetaken = testCountSort(arr5, ARRSIZE5);
    outputFile << timetaken << "\n";

    //testing Radix Sort
    outputFile << "RadixSort,";
    timetaken = testRadixSort(arr1, ARRSIZE1);
    outputFile << timetaken << ",";
    timetaken = testRadixSort(arr2, ARRSIZE2);
    outputFile << timetaken << ",";
    timetaken = testRadixSort(arr3, ARRSIZE3);
    outputFile << timetaken << ",";
    timetaken = testRadixSort(arr4, ARRSIZE4);
    outputFile << timetaken << ",";
    timetaken = testRadixSort(arr5, ARRSIZE5);
    outputFile << timetaken << "\n";

    cout << "The data has been saved to a file";
    //close file
    outputFile.close();
    

    return 0;
}

