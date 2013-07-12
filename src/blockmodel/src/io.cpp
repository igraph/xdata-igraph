/* vim:set ts=4 sw=4 sts=4 et: */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <block/blockmodel.h>
#include <block/io.hpp>
#include <block/util.hpp>
#include <igraph/cpp/graph.h>

using namespace igraph;
using namespace std;

/***************************************************************************/

typedef enum { SECTION_UNKNOWN, SECTION_INFO, SECTION_TYPES, SECTION_PROBABILITIES,
               SECTION_STICKINESSES, SECTION_RATES
} SectionType;

template<>
void PlainTextReader<UndirectedBlockmodel>::read(
        UndirectedBlockmodel* pModel, istream& is) {
    SectionType section = SECTION_UNKNOWN;
    map<string, int> headingMap;
    bool lookingForHeading = true;
    string str;
    stringstream iss;

    Vector types;
    Vector probabilities;
    m_filename.clear();

    headingMap["INFO"] = SECTION_INFO;
    headingMap["PROBABILITIES"] = SECTION_PROBABILITIES;
    headingMap["TYPES"] = SECTION_TYPES;

    if (is.fail())
        throw runtime_error("error while reading stream");

    while (!is.eof()) {
        getline(is, str);
        if (is.fail() && !is.eof())
            throw runtime_error("error while reading stream");

        /* Is the line empty? If so, a new heading will follow */
        if (str.size() == 0) {
            lookingForHeading = true;
            continue;
        }

        /* Are we looking for a heading? If so, store it */
        if (lookingForHeading) {
            map<string, int>::iterator it;
            it = headingMap.find(str);
            if (it != headingMap.end())
                section = static_cast<SectionType>(it->second);
            else
                section = SECTION_UNKNOWN;
            lookingForHeading = false;
            continue;
        }

        /* If we are here, we have the heading and the line is not empty */
        if (section == SECTION_INFO) {
            /* Store the original filename (if any) */
            string name;
            iss.str(str);
            iss >> name;
            if (name == "filename") {
                iss >> m_filename;
            }
        } else if (section == SECTION_TYPES) {
            /* Store the next type */
            int type;
            iss.str(str);
            iss >> type;
            if (iss.bad())
                throw runtime_error("error while parsing type vector");
            types.push_back(type);
        } else if (section == SECTION_PROBABILITIES) {
            /* Store the next row in the probability matrix */
            double p;
            iss.str(str);
            while (!iss.eof()) {
                iss >> p;
                if (iss.bad())
                    throw runtime_error("error while parsing probability matrix");
                probabilities.push_back(p);
            }
        }
        iss.clear();
    }

    // Find how many types we have
    size_t n = sqrt(probabilities.size());
    if (n * n != probabilities.size())
        throw runtime_error("probability matrix must be square");

    Matrix probMat(n, n);
    copy(probabilities.begin(), probabilities.end(), probMat.begin());

    pModel->setNumTypes(n);
    pModel->setTypes(types);
    pModel->setProbabilities(probMat);
}

template<>
void PlainTextReader<DegreeCorrectedUndirectedBlockmodel>::read(
        DegreeCorrectedUndirectedBlockmodel* pModel, istream& is) {
    SectionType section = SECTION_UNKNOWN;
    map<string, int> headingMap;
    bool lookingForHeading = true;
    string str;
    stringstream iss;

    Vector types;
    Vector thetas;
    Vector rates;
    m_filename.clear();

    headingMap["INFO"] = SECTION_INFO;
    headingMap["TYPES"] = SECTION_TYPES;
    headingMap["STICKINESSES"] = SECTION_STICKINESSES;
    headingMap["RATES"] = SECTION_RATES;

    if (is.fail())
        throw runtime_error("error while reading stream");

    while (!is.eof()) {
        getline(is, str);
        if (is.fail() && !is.eof())
            throw runtime_error("error while reading stream");

        /* Is the line empty? If so, a new heading will follow */
        if (str.size() == 0) {
            lookingForHeading = true;
            continue;
        }

        /* Are we looking for a heading? If so, store it */
        if (lookingForHeading) {
            map<string, int>::iterator it;
            it = headingMap.find(str);
            if (it != headingMap.end())
                section = static_cast<SectionType>(it->second);
            else
                section = SECTION_UNKNOWN;
            lookingForHeading = false;
            continue;
        }

        /* If we are here, we have the heading and the line is not empty */
        if (section == SECTION_INFO) {
            /* Store the original filename (if any) */
            string name;
            iss.str(str);
            iss >> name;
            if (name == "filename") {
                iss >> m_filename;
            }
        } else if (section == SECTION_TYPES) {
            /* Store the next type */
            int type;
            iss.str(str);
            iss >> type;
            if (iss.bad())
                throw runtime_error("error while parsing type vector");
            types.push_back(type);
        } else if (section == SECTION_STICKINESSES) {
            /* Store the next degree */
            double theta;
            iss.str(str);
            iss >> theta;
            if (iss.bad())
                throw runtime_error("error while parsing stickiness vector");
            thetas.push_back(theta);
        } else if (section == SECTION_RATES) {
            /* Store the next row in the rate matrix */
            double p;
            iss.str(str);
            while (!iss.eof()) {
                iss >> p;
                if (iss.bad())
                    throw runtime_error("error while parsing probability matrix");
                rates.push_back(p);
            }
        }
        iss.clear();
    }

    // Find how many types we have
    size_t n = sqrt(rates.size());
    if (n * n != rates.size())
        throw runtime_error("rate matrix must be square");

    Matrix rateMatrix(n, n);
    copy(rates.begin(), rates.end(), rateMatrix.begin());

    pModel->setNumTypes(n);
    pModel->setTypes(types);
    pModel->setStickinesses(thetas);
    pModel->setRates(rateMatrix);
}

void PlainTextReader<Blockmodel>::read(Blockmodel* model, istream& is) {
    if (typeid(*model) == typeid(UndirectedBlockmodel)) {
        PlainTextReader<UndirectedBlockmodel> reader;
        reader.read(static_cast<UndirectedBlockmodel*>(model), is);
        m_filename = reader.getFilename();
    } else if (typeid(*model) == typeid(DegreeCorrectedUndirectedBlockmodel)) {
        PlainTextReader<DegreeCorrectedUndirectedBlockmodel> reader;
        reader.read(static_cast<DegreeCorrectedUndirectedBlockmodel*>(model), is);
        m_filename = reader.getFilename();
    } else {
        throw std::runtime_error("reader does not know the given blockmodel");
    }
}

/***************************************************************************/

template<>
void PlainTextWriter<UndirectedBlockmodel>::write(
        const UndirectedBlockmodel* pModel, ostream& os) {
    const UndirectedBlockmodel& model = *pModel;
    const Graph* pGraph = model.getGraph();
    int n = pGraph->vcount();
    int k = model.getNumTypes();
    time_t now = time(0);

    os << "INFO\n";
    os << "date\t" << ctime(&now);
    if (pGraph->hasAttribute("filename"))
        os << "filename\t" << pGraph->getAttribute("filename").as<std::string>() << '\n';
    os << "num_vertices\t" << n << '\n';
    os << "num_types\t" << k << '\n';
    os << "log_likelihood\t" << model.getLogLikelihood() << '\n';
    os << "aic\t" << 2 * (model.getNumParameters() - model.getLogLikelihood()) << '\n';
    os << '\n';

    os << "TYPES\n";
    for (int i = 0; i < n; i++)
        os << i << '\t' << model.getType(i) << '\n';
    os << '\n';

    os << "PROBABILITIES\n";
    for (int i = 0; i < k; i++) {
        os << model.getProbability(i, 0);
        for (int j = 1; j < k; j++)
            os << '\t' << model.getProbability(i, j);
        os << '\n';
    }
}

template<>
void PlainTextWriter<DegreeCorrectedUndirectedBlockmodel>::write(
        const DegreeCorrectedUndirectedBlockmodel* pModel, ostream& os) {
    const DegreeCorrectedUndirectedBlockmodel& model = *pModel;
    const Graph* pGraph = model.getGraph();
    int n = pGraph->vcount();
    int k = model.getNumTypes();
    time_t now = time(0);

    os << "INFO\n";
    os << "date\t" << ctime(&now);
    if (pGraph->hasAttribute("filename"))
        os << "filename\t" << pGraph->getAttribute("filename").as<std::string>() << '\n';
    os << "num_vertices\t" << n << '\n';
    os << "num_types\t" << k << '\n';
    os << "log_likelihood\t" << model.getLogLikelihood() << '\n';
    os << "aic\t" << 2 * (model.getNumParameters() - model.getLogLikelihood()) << '\n';
    os << '\n';

    os << "TYPES\n";
    for (int i = 0; i < n; i++)
        os << model.getType(i) << '\n';
    os << '\n';

    os << "STICKINESSES\n";
    Vector thetas = model.getStickinesses();
    for (int i = 0; i < n; i++)
        os << thetas[i] << '\n';
    os << '\n';

    os << "RATES\n";
    for (int i = 0; i < k; i++) {
        os << model.getRate(i, 0);
        for (int j = 1; j < k; j++)
            os << '\t' << model.getRate(i, j);
        os << '\n';
    }
}

template<>
void PlainTextWriter<Blockmodel>::write(const Blockmodel* model, ostream& os) {
    if (typeid(*model) == typeid(UndirectedBlockmodel)) {
        PlainTextWriter<UndirectedBlockmodel> writer;
        writer.write(static_cast<const UndirectedBlockmodel*>(model), os);
    } else if (typeid(*model) == typeid(DegreeCorrectedUndirectedBlockmodel)) {
        PlainTextWriter<DegreeCorrectedUndirectedBlockmodel> writer;
        writer.write(static_cast<const DegreeCorrectedUndirectedBlockmodel*>(model), os);
    } else {
        throw std::runtime_error("writer does not know the given blockmodel");
    }
}

/***************************************************************************/

/// Converts a string to its JSON representation
static std::string toJSON(const std::string& str) {
    std::string out("\"");

    for (std::string::const_iterator it = str.begin(); it != str.end(); it++) {
        switch (*it) {
            case '"':
                out.append("\\\""); break;
            case '\\':
                out.append("\\\\"); break;
            case '\b':
                out.append("\\b"); break;
            case '\f':
                out.append("\\f"); break;
            case '\n':
                out.append("\\n"); break;
            case '\r':
                out.append("\\r"); break;
            case '\t':
                out.append("\\t"); break;
            case '/':
                out.append("\\/"); break;
            default:
                if (*it >= 0 && *it <= 31) {
                    char buf[7];
                    sprintf(buf, "\\u00%x", int(*it));
                    out.append(buf);
                } else out.append(1, *it);
        }
    }

    out.append("\"");

    return out;
}

template<>
void JSONWriter<UndirectedBlockmodel>::write(
        const UndirectedBlockmodel* pModel, ostream& os) {
    const UndirectedBlockmodel& model = *pModel;
    const Graph* pGraph = model.getGraph();
    int n = pGraph->vcount();
    int k = model.getNumTypes();
    time_t timestamp = time(0);
    char* formatted_date = ctime(&timestamp);

    // Trim the newline from the formatted date
    formatted_date[strlen(formatted_date)-1] = 0;

    os << "{\n"
       << "    \"info\": {\n"
       << "        \"date\": " << toJSON(formatted_date) << ",\n";

    if (pGraph->hasAttribute("filename"))
        os << "        \"filename\": "
           << toJSON(pGraph->getAttribute("filename").as<std::string>()) << ",\n";

    os << "        \"timestamp\": " << timestamp << ",\n"
       << "        \"num_vertices\": " << n << ",\n"
       << "        \"num_types\": " << k << ",\n"
       << "        \"log_likelihood\": " << model.getLogLikelihood() << ",\n"
       << "        \"aic\": " << aic(model) << "\n"
       << "    },\n"
       << "    \"parameters\": {\n"
       << "        \"types\": [";
    for (int i = 0; i < n-1; i++) {
        if (i % 15 == 0)
            os << "\n                  ";
        os << model.getType(i) << ", ";
    }
    if (n > 0)
        os << model.getType(n-1);
    os << "\n        ],\n";
    os << "        \"p\": [";
    if (k > 0) {
        for (int i = 0; i < k; i++) {
            os << "\n              [";
            for (int j = 0; j < k-1; j++) {
                os << model.getProbability(i, j) << ", ";
            }
            os << model.getProbability(i, k-1) << "]";
            if (i < k-1)
                os << ", ";
        }
    }
    os << "\n        ]\n"
       << "    }\n"
       << "}\n";
}

template<>
void JSONWriter<DegreeCorrectedUndirectedBlockmodel>::write(
        const DegreeCorrectedUndirectedBlockmodel* pModel, ostream& os) {
    const DegreeCorrectedUndirectedBlockmodel& model = *pModel;
    const Graph* pGraph = model.getGraph();
    int n = pGraph->vcount();
    int k = model.getNumTypes();
    time_t timestamp = time(0);
    char* formatted_date = ctime(&timestamp);

    // Trim the newline from the formatted date
    formatted_date[strlen(formatted_date)-1] = 0;

    os << "{\n"
       << "    \"info\": {\n"
       << "        \"date\": " << toJSON(formatted_date) << ",\n";

    if (pGraph->hasAttribute("filename"))
        os << "        \"filename\": "
           << toJSON(pGraph->getAttribute("filename").as<std::string>()) << ",\n";

    os << "        \"timestamp\": " << timestamp << ",\n"
       << "        \"num_vertices\": " << n << ",\n"
       << "        \"num_types\": " << k << ",\n"
       << "        \"log_likelihood\": " << model.getLogLikelihood() << ",\n"
       << "        \"aic\": " << aic(model) << "\n"
       << "    },\n"
       << "    \"parameters\": {\n"
       << "        \"types\": [";
    for (int i = 0; i < n-1; i++) {
        if (i % 15 == 0)
            os << "\n                  ";
        os << model.getType(i) << ", ";
    }
    if (n > 0)
        os << model.getType(n-1);
    os << "\n        ],\n"
       << "        \"stickinesses\": [";

    Vector stickinesses(model.getStickinesses());
    for (int i = 0; i < n-1; i++) {
        if (i % 15 == 0)
            os << "\n                  ";
        os << stickinesses[i] << ", ";
    }
    if (n > 0)
        os << stickinesses[n-1];
    os << "\n        ],\n";
    os << "        \"rates\": [";
    if (k > 0) {
        for (int i = 0; i < k; i++) {
            os << "\n              [";
            for (int j = 0; j < k-1; j++) {
                os << model.getRate(i, j) << ", ";
            }
            os << model.getRate(i, k-1) << "]";
            if (i < k-1)
                os << ", ";
        }
    }
    os << "\n        ]\n"
       << "    }\n"
       << "}\n";
}

template<>
void JSONWriter<Blockmodel>::write(const Blockmodel* model, ostream& os) {
    if (typeid(*model) == typeid(UndirectedBlockmodel)) {
        JSONWriter<UndirectedBlockmodel> writer;
        writer.write(static_cast<const UndirectedBlockmodel*>(model), os);
    } else if (typeid(*model) == typeid(DegreeCorrectedUndirectedBlockmodel)) {
        JSONWriter<DegreeCorrectedUndirectedBlockmodel> writer;
        writer.write(static_cast<const DegreeCorrectedUndirectedBlockmodel*>(model), os);
    } else {
        throw std::runtime_error("writer does not know the given blockmodel");
    }
}


