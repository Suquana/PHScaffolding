#include <iostream>   
#include <fstream>   
#include <string>    
#include <map>    
#include <vector>  
using namespace std;
struct readinfo {
	string reference;
	long long readstart;
	long long readend;
	char dir;
	long long referencestart;
	long long referenceend;
	long long length;
	readinfo(string rf, long long rds, long long rde, char d, long long rfs, long long rfe, long long len) :reference(rf), readstart(rds), readend(rde), dir(d), referencestart(rfs), referenceend(rfe), length(len) {}
};
int main(int argc, char* argv[]) {
	cerr << "-----------start smap----------------" << endl;
	map<string, vector<readinfo>> h_graph;
	if (argc != 3) {
		cerr << "用法: " << argv[0] << " <输入文件路径> <输出文件路径>" << endl;
		return 1;
	}
	string inpath = argv[1];
	string outpath = argv[2];
	ifstream inputfile(inpath);
	if (!inputfile.is_open()) {
		cerr << "无法打开文件: " << inpath << endl;
		return 1; // 退出程序表示错误
	}
	ofstream outputfile(outpath, ios::app);
	if (!outputfile.is_open()) {
		cerr << "无法打开文件: " << outpath << endl;
		return 1; // 退出程序表示错误
	}
	string line;
	long long ans = 0;
	cerr << "开始建图" << endl;
	while (getline(inputfile, line)) {
		ans++;
		if (ans % 1000 == 0) cerr << "已经读取了" << ans << "行数据" << endl;
		string s;
		string rf;
		string rd;
		long long rds, rde, rfs, rfe, len;
		char d;
		int ff = 0;
		bool fff = 0;
		long long num = 0;
		bool f = 0;
		for (long long i = 0; i < line.size(); i++) {
			if (line[i] == '	') {
				num++;
				if (num == 1) {
					auto it = h_graph.find(s);
					if (it != h_graph.end()) f = 1;
				}
				switch (num) {
				case 1:
					rd = s;
					break;
				case 2:
					break;
				case 3:
					rds = stoi(s);
					break;
				case 4:
					rde = stoi(s);
					break;
				case 5:
					d = s[0];
					break;
				case 6:
					rf = s;
					break;
				case 7:
					break;
				case 8:
					rfs = stoi(s);
					break;
				case 9:
					rfe = stoi(s);
					break;
				case 10:
					break;
				case 11:
					break;
				case 12:
					ff = stoi(s);
					break;
				}
				if (num == 12) break;
				s.clear();
			}
			else s += line[i];
			if (num == 12) break;
		}
		if (ff >= 50) {
			len = rde - rds;
			if (f == 0) {
				h_graph.insert({ rd , { readinfo(rf,rds,rde,d,rfs,rfe,len) } });
				//if (rf == "chr_8_contig_7457") cout << 1 << endl;
			}
			else {
				//if (rf == "chr_8_contig_7457") cout << 1 << endl;
				for (auto& ri : h_graph[rd]) {
					if (ri.reference == rf) {
						ri.length += len;
						fff = 1;
						break;
					}
				}
				if (fff == 0) {
					h_graph[rd].push_back(readinfo(rf, rds, rde, d, rfs, rfe, len));
				}
			}
		}
	}
	cerr << "建图完成" << endl;
	cerr << "开始写入" << endl;
	for (const auto& ele : h_graph) {
		string read = ele.first;
		outputfile << read << ':' << endl;
		//cout << read << endl;
		const vector<readinfo>& rinfo = ele.second;
		for (const readinfo& info : rinfo) {
			//cout << info.length << endl;
			//if (read == "000e5f4e-dc26-4838-b4a0-6badd612a002" && info.reference == "chr_4_contig_3490") cout << 1 << endl;
			if (info.length >= 100) {
				outputfile << info.reference << ',' << info.readstart << ',' << info.readend << ',' << info.dir << ',' << info.referencestart << ',' << info.referenceend << ',' << info.length << endl;
			} 
		}
	}
	cerr << "写入完成" << endl;
	inputfile.close();
	outputfile.close();
	return 0;
}