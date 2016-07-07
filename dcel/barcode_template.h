#ifndef __BARCODE_TEMPLATE_H__
#define __BARCODE_TEMPLATE_H__

#include <set>

struct BarTemplate
{
    unsigned begin;                   //index of xiMatrixEntry of the equivalence class corresponding to the beginning of this bar
    unsigned end;                     //index of xiMatrixEntry of the equivalence class corresponding to the end of this bar
    mutable unsigned multiplicity;    //maybe this is bad style, but multiplicity is not involved in comparisons

    BarTemplate(unsigned a, unsigned b);
    BarTemplate(unsigned a, unsigned b, unsigned m);
    BarTemplate(const BarTemplate& other);
    BarTemplate(); // for serialization

    bool operator<(const BarTemplate other) const;

    template <class Archive>
    void cerealize(Archive & ar) {
        ar(begin, end, multiplicity);
    }
};


class BarcodeTemplate
{
    public:
        BarcodeTemplate();

        void add_bar(unsigned a, unsigned b);               //adds a bar to the barcode template (updating multiplicity, if necessary)
        void add_bar(unsigned a, unsigned b, unsigned m);   //adds a bar with multiplicity to the barcode template

        std::set<BarTemplate>::iterator begin();    //returns an iterator to the first bar in the barcode
        std::set<BarTemplate>::iterator end();      //returns an iterator to the past-the-end element of the barcode
        bool is_empty();                            //returns true iff this barcode has no bars

        void print();   //for testing only

        template <class Archive>
                void cerealize(Archive &ar) {
            ar(bars);
        }
    private:
        std::set<BarTemplate> bars;
};

#endif // __BARCODE_TEMPLATE_H__
