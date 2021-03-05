// Sequences-negative.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <set>
#include <map>

int length = 0;
int candidatecurl = 2;
int candidateperiod = 1;
std::vector<int> Tail = {};
std::vector<int> Periods = {};
std::vector<int> Generator = {};
std::vector<int> Max_tail_lengths = {};
std::vector<std::vector<int>> Generators_memory = {};
std::vector<std::vector<int>> Best_generators = {};
std::set<int> Change_indices = {};
std::map<int, std::set<int>> Dict = {};
std::map<int, std::map<int, std::set<int>>> Dicts_memory = {};

std::vector<int> seq_new = {};
std::map<int, std::set<int>> Dict_new = {};

void krul(std::vector<int>* seq, int* curl, int* period) {           // curl = 1, period = 0
    int l = seq->size();
    for (int i = 1; i <= (l / 2); ++i) {
        int j = i;
        while ((*seq)[l - j - 1] == (*seq)[l - j - 1 + i]) {
            ++j;
            if (j >= l) 
                break;
        }
        int candidate = j / i;
        if (candidate > *curl) {
            *curl = candidate;
            *period = i;
        }
    }
}

void tail_with_periods(std::vector<int> seq, std::vector<int>* tail, std::vector<int>* periods) {               // tail = {}, periods = {}
    int curl = 1, period = 0;
    std::vector<int> temp = seq;
    krul(&seq, &curl, &period);
    while (curl > 1) {
        tail->push_back(curl);
        temp.push_back(curl);           // easy optimization
        periods->push_back(period);
        curl = 1, period = 0;
        krul(&temp, &curl, &period);
    }
}

void tail_with_periods_part(std::vector<int> seq, std::vector<int>* tail, std::vector<int>* periods, int i) {   // tail = {}, periods = {}
    int curl = 1, period = 0;
    std::vector<int> temp = seq;
    krul(&seq, &curl, &period);
    while (curl > 1 and tail->size() < i) {
        tail->push_back(curl);
        temp.push_back(curl);           // easy optimization
        periods->push_back(period);
        curl = 1, period = 0;
        krul(&temp, &curl, &period);
    }
}

bool check_period_size() {
    return (candidatecurl * candidateperiod) > (length + Periods.size());
}

bool check_candidatecurl_size() {
    if (Tail.size())
        return candidatecurl > Tail.back() + 1;
    else
        return candidatecurl > length;
}

void up() {
    ++candidateperiod;
    while (check_period_size()) {
        ++candidatecurl;
        candidateperiod = 1;
        if (check_candidatecurl_size) {
            auto index = std::find(Change_indices.begin(), Change_indices.end(), Periods.size());
            if (index != Change_indices.end()) {
                Change_indices.erase(index);
            }
            if (Tail.size() == 0) {
                candidatecurl = 0;
                break;
            }
            index = std::find(Change_indices.begin(), Change_indices.end(), Periods.size() - 1);
            if (index == Change_indices.end()) {
                Change_indices.insert(Periods.size() - 1);
                Dict[Tail.back()].erase(length + Tail.size() - 1);
                candidatecurl = Tail.back() + 1;
                candidateperiod = 1;
                Tail.pop_back();
                Periods.pop_back();
            }
            else {
                candidatecurl = Tail.back();
                candidateperiod = Periods.back() + 1;
                Tail.pop_back();
                Periods.pop_back();
                Generator = Generators_memory[Periods.size()];
                Dict = Dicts_memory[Periods.size()];
                auto it = Generators_memory.begin();
                it->erase(it->begin() + Periods.size());
                Dicts_memory.erase(Periods.size());
            }
        }
    }
}

int real_generator_length() {
    int i = 0;
    while (Dict[Generator[i]].size() == 1) {
        ++i;
        if (i == length)
            break;
    }
    return length - i;
}

bool check_positive(int len) {
    for (int i : std::vector<int>(Generator.begin() + length - len, Generator.end())) {
        if (i < 1)
            return false;
    }
    return true;
}

void append() {
    Generators_memory[Periods.size()] = Generator;
    Dicts_memory[Periods.size()] = Dict;
    Generator = std::vector<int>(seq_new.begin(), seq_new.begin() + length);
    Dict = Dict_new;
    Dict[candidatecurl].insert(length + Tail.size());
    Tail.push_back(candidatecurl);
    Periods.push_back(candidateperiod);

    int curl = 1, period = 0;
    std::vector<int> temp = Generator;
    temp.push_back(candidatecurl);
    while (true) {
        curl = 1, period = 0;
        krul(&temp, &curl, &period);
        if (curl == 1)
            break;
        Tail.push_back(curl);
        temp.push_back(curl);
        Periods.push_back(period);
        Dict[curl].insert(length + Tail.size() - 1);
    }
    candidatecurl = 2;
    candidateperiod = 1;
    Change_indices.insert(Tail.size());
    int len = real_generator_length();
    if (Max_tail_lengths.back() == Tail.size()) {
        std::vector<int> temp = Best_generators.back();
        temp.insert(temp.end(), Generator.begin(), Generator.end() - len);
        Best_generators.back() = temp;
    }
    if (Max_tail_lengths.back() < Tail.size()) {
        Max_tail_lengths.back() = Tail.size();
        Best_generators.back() = std::vector<int>(Generator.begin(), Generator.end() - len);
    }
}

bool test_1() {
    seq_new = Generator;
    seq_new.insert(seq_new.end(), Tail.begin(), Tail.end());
    Dict_new = Dict;
    int l = seq_new.size();
    for (int i = 0; i < (candidatecurl - 1) * candidateperiod; ++i) {
        int a = seq_new[l - 1 - i];
        int b = seq_new[l - 1 - i - candidateperiod];
        if (a != b and a > 0 and b > 0)
            return false;
        if (a > b) {
            for (int j = 0; j < l; ++j) {
                if (seq_new[j] == b)
                    seq_new[j] = a;
            }
            for (int x : Dict_new[b])
                Dict_new[a].insert(x);
            Dict_new.erase(b);
        }
        else if (a < b) {
            for (int j = 0; j < l; ++j) {
                if (seq_new[j] == a)
                    seq_new[j] = b;
            }
            for (int x : Dict_new[a])
                Dict_new[b].insert(x);
            Dict_new.erase(a);
        }
    }
    return true;
}

bool test_2() {
    int l = seq_new.size();
    std::vector<int> temp_seq = seq_new;
    temp_seq.push_back(candidatecurl);
    std::vector<int> temp_period = Periods;
    temp_period.push_back(candidateperiod);
    int curl = 1, period = 0;
    for (int i = 0; i <= length + l; ++i) {
        std::vector<int> temp = std::vector<int>(seq_new.begin(), seq_new.begin() + length + i);
        curl = 1, period = 0;
        krul(&temp, &curl, &period);
        if (curl != temp_seq[length + i] or period != temp_period[i])
            return false;
    }
    return true;
}

bool check_if_period_works() {
    if (test_1()) {
        if (test_2())
            return true;
    }
    return false;
}

void backtracking_step() {
    if (check_if_period_works())
        append();
    else
        up();
}

void backtracking(int k1, int p1, int k2, int p2) {
    Change_indices.insert(0);
    for (int i = 0; i < length; ++i) {
        Generator.push_back(-i);
        Dict[-i].insert(i);
        Max_tail_lengths.push_back(0);
        Best_generators.push_back({ 0 });
    }
    std::vector<int> seq = Generator;
    seq.insert(seq.end(), Tail.begin(), Tail.end());

    while (candidatecurl) {
        if (Tail.size() == 0 and candidatecurl == k2 and candidateperiod == p2)
            break;
        backtracking_step();
    }
    int record = 0;
    for (int i = 0; i < length; ++i) {
        if (Max_tail_lengths[i] > record) {
            record = Max_tail_lengths[i];
            std::cout << i + 1 << ": " << record << ", [";
            for (int x : Best_generators[i])
                std::cout << x << ", ";
            std::cout << std::endl;
        }
    }
}

int main()
{
    int k1, k2, p1, p2;
    std::cout << "Length: ";
    std::cin >> length;
    std::cin >> k1;
    std::cin >> p1;
    std::cin >> k2;
    std::cin >> p2;
    backtracking(k1, p1, k2, p2);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
