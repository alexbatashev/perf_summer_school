#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

class EventID {
public:
  EventID() : mId(-1) {}

private:
  friend class EventDispatcher;

  EventID(unsigned char id) : mId(id) {}

  unsigned char mId;
};

class EventDispatcher {
public:
  EventDispatcher() = default;

  EventID acquire() {
    size_t idx = __builtin_ffs(mFreeList);
    mFreeList |= 1 << idx;
    return EventID(idx);
    // for (std::size_t i = 0; i < mFreeList.size(); i++) {
    //   if (mFreeList[i] == false) {
    //     mFreeList[i] = true;
    //     return EventID(i);
    //   }
    // }
    //
    // return EventID();
  }

  void release(EventID id) {
    if (id.mId < 64) {
      mFreeList &= ~(1 << id.mId);
    }
    // mFreeList[id.mId] = false;
  }

  bool hasFreeSlots() {
    return mFreeList != std::numeric_limits<uint64_t>::max();
    // for (bool b : mFreeList) {
    //   if (b == false) return true;
    // }
    //
    // return false;
  }
private:
  // std::array<bool, 64> mFreeList = {false};
  uint64_t mFreeList = 0;
};

